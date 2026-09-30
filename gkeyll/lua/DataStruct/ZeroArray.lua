-- Gkyl ------------------------------------------------------------------------
--
-- Lua interface to gkylzero's gkyl_array structure
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local Alloc = require "Lib.Alloc"
local cuda
require "Lib.ZeroUtil"

_M = {}

ffi.cdef [[
/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
  enum gkyl_elem_type type; // type of data stored in array
  size_t elemsz, ncomp; // size of elements, number of 'components'
  size_t size; // number of indices

  size_t esznc; // elemsz*ncomp
  void *_data; // pointer to data
  uint32_t flags;  
  struct gkyl_ref_count ref_count;

  int nthreads, nblocks; // threads per block, number of blocks
  struct gkyl_array *on_dev; // pointer to itself or device data
};

// To return diff of two arrays
struct gkyl_array_diff {
  bool is_compatible; // are arrays compatible

  // the following make sense only if is_compatible = true
  double max_abs_diff; // maximum absolute difference
  double min_abs_diff; // minmum absolute difference
  double max_rel_diff; // maximum relative difference
  double min_rel_diff; // minmum relative difference
};  

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with data on NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the
 * host.  However, the on_dev member for this cal is set to a device
 * clone of the host struct, and is what should be used to pass to
 * CUDA kernels which require the entire array struct on device.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with host-pinned data for use with NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Returns true if array lives on NV-GPU.
 *
 * @param arr Array to check
 * @return true of array lives on NV-GPU, false otherwise
 */
bool gkyl_array_is_cu_dev(const struct gkyl_array *arr);

enum gkyl_array_op { GKYL_MIN, GKYL_MAX, GKYL_SUM, GKYL_SUM_ABS };


/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Copy into array using async methods for cuda arrays: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy_async(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear(struct gkyl_array *out, double val);

/**
 * Compute out = out + a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Compute out = out + a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = out[coff]+ a*inp if
 * out->ncomp > inp->ncomp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_offset(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_set(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Set out = a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = a*inp if
 * out->ncomp > inp->ncomp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_set_offset(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff);

/**
 * Scale out = a*out. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale
 * @return out array
 */
struct gkyl_array* gkyl_array_scale(struct gkyl_array *out, double a);

struct gkyl_array* gkyl_array_scale_by_cell(struct gkyl_array *out, const struct gkyl_array *a);

/**
 * Shift the k-th coefficient in every cell, out_k = a+out_k. Returns out.
 *
 * @param out Output array.
 * @param a Factor to shift k-th coefficient by.
 * @param k Coefficient to be shifted.
 * @return out array.
 */
struct gkyl_array* gkyl_array_shiftc(struct gkyl_array *out, double a, unsigned k);

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear_range(struct gkyl_array *out, double val,
  const struct gkyl_range *range);

/**
 * Compute out = out + a*inp over a range of indices.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param range Range specifying region to accumulate
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Compute out = out + a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = out[coff]+ a*inp if
 * out->ncomp > inp->ncomp, over a range of indices. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Set out = a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = a*inp if
 * out->ncomp > inp->ncomp, over a range of indices. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range);

/**
 * Scale out = a*ut. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale by
 * @return out array
 * @param range Range specifying region to scale
 */
struct gkyl_array* gkyl_array_scale_range(struct gkyl_array *out,
  double a, const struct gkyl_range *range);

/**
 * Shift the k-th coefficient in every cell, out_k = a+out_k within
 * a given range. Returns out.
 *
 * @param out Output array.
 * @param a Factor to shift k-th coefficient by.
 * @param k Coefficient to be shifted.
 * @param range Range to shift coefficient k in.
 * @return out array.
 */
struct gkyl_array* gkyl_array_shiftc_range(struct gkyl_array *out, double a,
  unsigned k, const struct gkyl_range *range);

/**
 * Shift the k-th coefficient in every cell, out_k = a+out_k within
 * a given range. Returns out.
 *
 * @param out Output array.
 * @param a Factor to shift k-th coefficient by.
 * @param k Coefficient to be shifted.
 * @param range Range to shift coefficient k in.
 * @return out array.
 */
struct gkyl_array* gkyl_array_shiftc_range(struct gkyl_array *out, double a,
  unsigned k, struct gkyl_range range);

/**
 * Copy out inp. Returns out.
 *
 * @param out Output array
 * @param inp Input array
 * @param range Range specifying region to copy
 * @return out array
 */
struct gkyl_array* gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Copy out inp over specified ranges. Returns out.
 * input and output ranges must have the same volume.
 *
 * @param out Output array
 * @param inp Input array
 * @param out_range Range specifying region to copy to in out array
 * @param inp_range Range specifying region to copy to from in inp array
 * @return out array
 */
struct gkyl_array* gkyl_array_copy_range_to_range(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range *out_range, struct gkyl_range *inp_range);

/**
 * Perform an "reduce" operation of data in the array.
 *
 * @param res On output, reduces values (ncomp size)
 * @param arr Array to perform reduction on.
 * @param op Reduction operators
 */
void gkyl_array_reduce(double *res, const struct gkyl_array *arr, enum gkyl_array_op op);

/**
 * Perform an "reduce" operation of data in the array. Data is reduced
 * component-wise.
 *
 * @param res On output, reduced values (ncomp size)
 * @param arr Array to perform reduction on.
 * @param op Reduction operators
 * @param range Range specifying region
 */
void gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);

/**
 * Copy region of array into a buffer. The buffer must be preallocated
 * and at least of size arr->size*arr->elemSz bytes.
 *
 * @param arr Array to copy from
 * @param data Output data buffer.
 * @param range Range specifying region to copy from
 */
void gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  struct gkyl_range *range);

/**
 * Copy buffer into region of array. The array must be preallocated.
 *
 * @param arr Array to copy into
 * @param data Input data buffer.
 * @param range Range specifying region to copy into
 */
void gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, struct gkyl_range *range);

/**
 * Release pointer to array
 *
 * @param arr Array to release.
 */
void gkyl_array_release(const struct gkyl_array* arr);

/**
 * Return difference between two arrays. Mostly useful for testing.
 *
 * @param arr1 First array to compare
 * @param arr2 Second array to compare
 * @param range Range to compare over
 * @return diff between arrays
 */
struct gkyl_array_diff gkyl_array_diff(const struct gkyl_array *arr1,
  const struct gkyl_array *arr2, const struct gkyl_range *range);

]]

local longSz = sizeof("long")

-- Various types for arrays of basic C-types
_M.int    = 'GKYL_INT'
_M.int64  = 'GKYL_INT64'
_M.float  = 'GKYL_FLOAT'
_M.double = 'GKYL_DOUBLE'
_M.long   = typeof('long')
_M.char   = typeof('char')

local function getArrayTypeCode(atype)
   if atype == typeof("int") then
      return 1
   elseif atype == typeof("int64") then
      return 2
   elseif atype == typeof("float") then
      return 3
   elseif atype == typeof("double") then
      return 4
   end
   return 42 -- user-defined type
end

local function getType(enum)
   if enum == 0 then
      return "int"
   elseif enum == 1 then
      return "int64"
   elseif enum == 2 then
      return "float"
   elseif enum == 3 then
      return "double"
   end
end

-- Array ctype
local ArrayCt = typeof("struct gkyl_array")

local array_fn = {
   copy = function (self, src)
      return ffiC.gkyl_array_copy(self, src)
   end,
   copyAsync = function (self, src)
      return ffiC.gkyl_array_copy_async(self, src)
   end,
   clone = function (self)
      return ffiC.gkyl_array_clone(self)
   end,
   aquire = function (self)
      return ffiC.gkyl_array_acquire(self)
   end,
   data = function (self)
      return self:fetch(0)
   end,
   fetch = function (self, loc)
      if loc == nil then loc = 0 end
      return ffi.cast(getType(self.type).."*", self._data) + loc*self.ncomp
   end,
   cfetch = function (self, loc)
      if loc == nil then loc = 0 end
      return ffi.cast("const "..getType(self.type).."*", self._data) + loc*self.ncomp
   end,
   get_size = function (self)
      return tonumber(self.size)
   end,
   get_ncomp = function (self)
      return tonumber(self.ncomp)
   end,
   get_elemsz = function (self)
      return tonumber(self.elemsz)
   end,
   clear = function (self, val)
      ffiC.gkyl_array_clear(self, val)
   end,
   clearRange = function (self, val, rng)
      ffiC.gkyl_array_clear_range(self, val, rng)
   end,
   set = function (self, val, fld)
      ffiC.gkyl_array_set(self, val, fld)
   end,
   setOffset = function (self, val, fld, off)
      ffiC.gkyl_array_set_offset(self, val, fld, off)
   end,
   accumulate = function (self, val, fld)
      ffiC.gkyl_array_accumulate(self, val, fld)
   end,
   accumulateOffset = function (self, val, fld, off)
      ffiC.gkyl_array_accumulate_offset(self, val, fld, off)
   end,
   setRange = function (self, val, fld, rng)
      ffiC.gkyl_array_set_range(self, val, fld, rng)
   end,
   setOffsetRange = function (self, val, fld, off, rng)
      ffiC.gkyl_array_set_offset_range(self, val, fld, off, rng)
   end,
   copyRangeToRange = function (self, infld, outrng, inrng)
      ffiC.gkyl_array_copy_range_to_range(self, infld, outrng, inrng)
   end,
   copyRange = function (self, infld, inrng)
      ffiC.gkyl_array_copy_range(self, infld, inrng)
   end,
   accumulateRange = function (self, val, fld, rng)
      ffiC.gkyl_array_accumulate_range(self, val, fld, rng)
   end,
   accumulateOffsetRange = function (self, val, fld, off, rng)
      ffiC.gkyl_array_accumulate_offset_range(self, val, fld, off, rng)
   end,
   scale = function (self, val)
      ffiC.gkyl_array_scale(self, val)
   end,
   scale_by_cell = function (self, val)
      ffiC.gkyl_array_scale_by_cell(self, val)
   end,
   shiftc = function (self, val, comp)
      ffiC.gkyl_array_shiftc(self, val, comp)
   end,
   shiftcRange = function (self, val, comp, rng)
      ffiC.gkyl_array_shiftc_range(self, val, comp, rng)
   end,
   reduceRange = function (self, out, op, rng)
      if op == "min" then 
         enum = 0 
      elseif op == "max" then
         enum = 1
      elseif op == "sum" then
         enum = 2
      elseif op == "sum_abs" then
         enum = 3
      end
      ffiC.gkyl_array_reduce_range(out, self, enum, rng)
   end,
   copy_to_buffer = function (self, data, rng)  
      ffiC.gkyl_array_copy_to_buffer(data, self, rng)
   end,
   copy_from_buffer = function (self, data, rng)  
      ffiC.gkyl_array_copy_from_buffer(self, data, rng)
   end,
   is_cu_dev = function (self)
      return ffiC.gkyl_array_is_cu_dev(self)
   end,
}

local array_mt = {
   __new = function (self, atype, ncomp, size, on_gpu)
      if on_gpu==1 then
         return ffi.gc(ffiC.gkyl_array_cu_dev_new(atype, ncomp, size),
                       ffiC.gkyl_array_release)
      elseif on_gpu==2 then
         return ffi.gc(ffiC.gkyl_array_cu_host_new(atype, ncomp, size),
                       ffiC.gkyl_array_release)
      else
         return ffi.gc(ffiC.gkyl_array_new(atype, ncomp, size),
                       ffiC.gkyl_array_release)
      end
   end,
   __index = array_fn,
}
local ArrayCtor = metatype(ArrayCt, array_mt)

-- Construct array of given shape and type
_M.Array = function (atype, ncomp, size, on_gpu)
   return ArrayCtor(atype, ncomp, size, on_gpu)
end

-- Compare arrays a1 nd a2 and return struct gkyl_array_diff
_M.arrayDiff = function(a1, a2, range)
   return ffi.C.gkyl_array_diff(a1, a2, range)
end

return _M

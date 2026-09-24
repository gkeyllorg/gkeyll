// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifdef GKYL_HAVE_MPI

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include <gkyl_mpi_funcs.h>

// Sizeof operators for various objects
GET_MPI_OBJ_SIZE(MPI_Status);
GET_MPI_OBJ_SIZE(MPI_Aint);
GET_MPI_OBJ_SIZE(MPI_Request);
GET_MPI_OBJ_PTR_SIZE(MPI_Status);
GET_MPI_OBJ_PTR_SIZE(MPI_Request);

// Pre-defined objects and constants
GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);
GET_MPI_OBJECT(Comm, MPI_COMM_NULL);
GET_MPI_OBJECT(Comm, MPI_COMM_SELF)
GET_MPI_OBJECT(Request, MPI_REQUEST_NULL);
GET_MPI_OBJECT_PTR(Status, MPI_STATUS_IGNORE);
GET_MPI_OBJECT(Info, MPI_INFO_NULL);
GET_INT_OBJECT(MPI_PROC_NULL);

// Datatypes
GET_MPI_OBJECT(Datatype, MPI_C_BOOL);
GET_MPI_OBJECT(Datatype, MPI_CHAR);
GET_MPI_OBJECT(Datatype, MPI_BYTE);
GET_MPI_OBJECT(Datatype, MPI_SHORT);
GET_MPI_OBJECT(Datatype, MPI_INT);
GET_MPI_OBJECT(Datatype, MPI_LONG);
GET_MPI_OBJECT(Datatype, MPI_FLOAT);
GET_MPI_OBJECT(Datatype, MPI_DOUBLE);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_CHAR);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_SHORT);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_LONG);
GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE);
GET_MPI_OBJECT(Datatype, MPI_LONG_LONG_INT);
GET_MPI_OBJECT(Datatype, MPI_FLOAT_INT);
GET_MPI_OBJECT(Datatype, MPI_LONG_INT);
GET_MPI_OBJECT(Datatype, MPI_DOUBLE_INT);
GET_MPI_OBJECT(Datatype, MPI_SHORT_INT);
GET_MPI_OBJECT(Datatype, MPI_2INT);
GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE_INT);
GET_MPI_OBJECT(Datatype, MPI_PACKED);

// Operators
GET_MPI_OBJECT(Op, MPI_MAX);
GET_MPI_OBJECT(Op, MPI_MIN);
GET_MPI_OBJECT(Op, MPI_SUM);
GET_MPI_OBJECT(Op, MPI_PROD);
GET_MPI_OBJECT(Op, MPI_LAND);
GET_MPI_OBJECT(Op, MPI_BAND);
GET_MPI_OBJECT(Op, MPI_LOR);
GET_MPI_OBJECT(Op, MPI_BOR);
GET_MPI_OBJECT(Op, MPI_LXOR);
GET_MPI_OBJECT(Op, MPI_BXOR);
GET_MPI_OBJECT(Op, MPI_MINLOC);
GET_MPI_OBJECT(Op, MPI_MAXLOC);

// Error codes
GET_INT_OBJECT(MPI_SUCCESS);

// Constants
GET_INT_OBJECT(MPI_COMM_TYPE_SHARED);
GET_INT_OBJECT(MPI_UNDEFINED);
GET_INT_OBJECT(MPI_ORDER_C);
GET_INT_OBJECT(MPI_ORDER_FORTRAN);
GET_INT_OBJECT_PTR(MPI_UNWEIGHTED);
GET_VOID_OBJECT_PTR(MPI_BOTTOM);

// Functions to allocate and free structs holding requests and statuses.
void gkyl_MPI_Request_alloc(gkyl_MPI_Request *rs, int num) {
  rs->req = (MPI_Request *) malloc(num*sizeof(MPI_Request));
}
void gkyl_MPI_Request_release(gkyl_MPI_Request *rs) {
  free(rs->req);
}
void gkyl_MPI_Status_alloc(gkyl_MPI_Status *ss, int num) {
  ss->stat = (MPI_Status *) malloc(num*sizeof(MPI_Status));
}
void gkyl_MPI_Status_release(gkyl_MPI_Status *ss) {
  free(ss->stat);
}
void gkyl_MPI_Request_Status_alloc(gkyl_MPI_Request_Status *rss, int num) {
  rss->req = (MPI_Request *) malloc(num*sizeof(MPI_Request));
  rss->stat = (MPI_Status *) malloc(num*sizeof(MPI_Status));
}
void gkyl_MPI_Request_Status_release(gkyl_MPI_Request_Status *rss) {
  free(rss->req);
  free(rss->stat);
}

// Functions to fetch members of status.
int gkyl_mpi_get_status_SOURCE(const MPI_Status* instat, int off) {
  return instat[off].MPI_SOURCE;
}
int gkyl_mpi_get_status_TAG(const MPI_Status* instat, int off) {
  return instat[off].MPI_TAG;
}
int gkyl_mpi_get_status_ERROR(const MPI_Status* instat, int off) {
  return instat[off].MPI_ERROR;
}

// Get count from a status (which may be one of several in an array of
// statuses).
int gkyl_mpi_get_status_count(const MPI_Status *instat, MPI_Datatype datatype, int *count, int off) {
  int err = MPI_Get_count(instat+off, datatype, count);
  return err;
}

#endif

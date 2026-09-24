/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_array_ops.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>

  struct gkyl_array_container {
    struct gkyl_array *arr;
  };
  
  struct gkyl_container_pack {
    struct gkyl_array_container *ac;
  };

  struct gkyl_array_bag {
    struct gkyl_array *arr;
    struct gkyl_array_bag *bag;
  };

  // Functions for test_array_container.
  void test_array_container_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_containers,
    struct gkyl_array_container *acs1, struct gkyl_array_container *acs2);
  void test_array_container_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_containers,
    struct gkyl_array_container *acs1, double a, struct gkyl_array_container *acs2);
  int test_array_container_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_containers,
    struct gkyl_array_container *acs1);

  // Functions for test_array_bag.
  void test_array_bag_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_bags,
    struct gkyl_array_bag *bag1, struct gkyl_array_bag *bag2);
  void test_array_bag_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_bags,
    struct gkyl_array_bag *bag1, double a, struct gkyl_array_bag *bag2);
  int test_array_bag_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_bags,
    struct gkyl_array_bag *bag1);
}

#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)
// Compute number of elements stored in array 'arr'
#define NELM(arr) (arr->size*arr->ncomp)
// Compute size of 'arr'
#define NSIZE(arr) (arr->size)
// Compute number of components stored in array 'arr'
#define NCOM(arr) (arr->ncomp)

//
// Functions for test_array_container.
//

__global__
void ker_cu_array_container_accumulate_dev_assign(int num_containers,
  struct gkyl_array_container *acs1, struct gkyl_array_container *acs2)
{
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k];
    double *arr1_d = (double*) arrc1->arr->data;
    for (unsigned long linc = START_ID; linc < NELM(arrc1->arr); linc += blockDim.x*gridDim.x) {
      arr1_d[linc] = k*100.0 + linc*1.0;
    }

    struct gkyl_array_container *arrc2 = &acs2[k];
    double *arr2_d  = (double*) arrc2->arr->data;
    for (unsigned long linc = START_ID; linc < NELM(arrc2->arr); linc += blockDim.x*gridDim.x) {
      arr2_d[linc] = k*200.0 + linc*2.0;
    }
  }
}

void
test_array_container_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1, struct gkyl_array_container *acs2)
{
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_container_accumulate_dev_assign<<<nblocks,nthreads>>>(num_containers, acs1, acs2);
}

__global__
void ker_cu_array_container_accumulate_dev_accumulate(int num_containers, 
  struct gkyl_array_container *acs1, double a, struct gkyl_array_container *acs2)
{
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
    double *arr1_d = (double*) arrc1->arr->data;
    double *arr2_d  = (double*) arrc2->arr->data;
    for (unsigned long linc = START_ID; linc < NELM(arrc1->arr); linc += blockDim.x*gridDim.x) {
      arr1_d[linc] += a*arr2_d[linc];
    }
  }
}

void
test_array_container_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1, double a, struct gkyl_array_container *acs2)
{
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_container_accumulate_dev_accumulate<<<nblocks,nthreads>>>(num_containers, acs1, a, acs2);
}

__global__
void ker_cu_array_container_accumulate_dev_check(int num_containers, 
  struct gkyl_array_container *acs1, int *nfail)
{
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k];
    double *arr1_d = (double*) arrc1->arr->data;
    for (unsigned long linc = START_ID; linc < NELM(arrc1->arr); linc += blockDim.x*gridDim.x) {
      GKYL_CU_CHECK( arr1_d[linc] == 2.0*(k*100.0+linc*1.0), nfail );
    }
  }
}

int
test_array_container_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));

  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_container_accumulate_dev_check<<<nblocks,nthreads>>>(num_containers, acs1, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

//
// Functions for test_array_bag.
//

__global__
void ker_cu_array_bag_accumulate_dev_assign(int num_arrays,
  struct gkyl_array_bag *bag1, struct gkyl_array_bag *bag2)
{
  for (int k=0; k<num_arrays; k++) {
    struct gkyl_array_bag *innerbag1 = &bag1[k];
    struct gkyl_array *arr1 = innerbag1->arr;
    double *arr1_d = (double *) arr1->data;
    for (unsigned long linc = START_ID; linc < NELM(arr1); linc += blockDim.x*gridDim.x) {
      arr1_d[linc] = k*100.0 + linc*1.0;
    }

    struct gkyl_array_bag *innerbag2 = &bag2[k];
    struct gkyl_array *arr2 = innerbag2->arr;
    double *arr2_d = (double *) arr2->data;
    for (unsigned long linc = START_ID; linc < NELM(arr2); linc += blockDim.x*gridDim.x) {
      arr2_d[linc] = k*200.0 + linc*2.0;
    }
  }
}

void
test_array_bag_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_arrays,
  struct gkyl_array_bag *innerbag1, struct gkyl_array_bag *innerbag2)
{
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_bag_accumulate_dev_assign<<<nblocks,nthreads>>>(num_arrays, innerbag1, innerbag2);
}

__global__
void ker_cu_array_bag_accumulate_dev_accumulate(int num_bags, 
  struct gkyl_array_bag *bag1, double a, struct gkyl_array_bag *bag2)
{
  for (int k=0; k<num_bags; k++) {
    struct gkyl_array_bag *innerbag1 = &bag1[k];
    struct gkyl_array *arr1 = innerbag1->arr;
    double *arr1_d = (double *) arr1->data;

    struct gkyl_array_bag *innerbag2 = &bag2[k];
    struct gkyl_array *arr2 = innerbag2->arr;
    double *arr2_d = (double *) arr2->data;

    for (unsigned long linc = START_ID; linc < NELM(arr1); linc += blockDim.x*gridDim.x) {
      arr1_d[linc] += a*arr2_d[linc];
    }
  }
}

void
test_array_bag_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_bags,
  struct gkyl_array_bag *innerbag1, double a, struct gkyl_array_bag *innerbag2)
{
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_bag_accumulate_dev_accumulate<<<nblocks,nthreads>>>(num_bags, innerbag1, a, innerbag2);
}

__global__
void ker_cu_array_bag_accumulate_dev_check(int num_bags, 
  struct gkyl_array_bag *bag1, int *nfail)
{
  for (int k=0; k<num_bags; k++) {
    struct gkyl_array_bag *innerbag1 = &bag1[k];
    struct gkyl_array *arr1 = innerbag1->arr;
    double *arr1_d = (double *) arr1->data;
    for (unsigned long linc = START_ID; linc < NELM(arr1); linc += blockDim.x*gridDim.x) {
      GKYL_CU_CHECK( arr1_d[linc] == 2.0*(k*100.0+linc*1.0), nfail );
    }
  }
}

int
test_array_bag_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_bags,
  struct gkyl_array_bag *innerbag1)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));

  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(arr_size*arr_ncomp, nthreads);
  ker_cu_array_bag_accumulate_dev_check<<<nblocks,nthreads>>>(num_bags, innerbag1, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

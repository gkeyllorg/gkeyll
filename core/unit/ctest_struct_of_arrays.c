#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

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

void test_array_container_accumulate_ho()
{
  // Test an approach to creating an array of arrays using the gkyl_array_container struct.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_containers = 2; // Number of array containers (i.e. no. of arrays).

  // Allocate objects.
  struct gkyl_array_container *acs1 = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
  struct gkyl_array_container *acs2 = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));

  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
    arrc1->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    arrc2->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
  }

  // Assign arrays.
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
    double *arr1_d = arrc1->arr->data;
    for (unsigned i=0; i<arrc1->arr->size; ++i) {
      arr1_d[i] = k*100.0 + i*1.0;
    }

    double *arr2_d  = arrc2->arr->data;
    for (unsigned i=0; i<arrc2->arr->size; ++i) {
      arr2_d[i] = k*200.0 + i*2.0;
    }
  }

  // Accumulate arrays.
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
    gkyl_array_accumulate(arrc1->arr, 0.5, arrc2->arr);
  }

  // Check results.
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k];
    double *arr1_d  = arrc1->arr->data;
    for (unsigned i=0; i<arrc1->arr->size; ++i)
      TEST_CHECK( gkyl_compare(arr1_d[i], 2.0*(k*100.0+i*1.0), 1e-14) );
  }

  // Free objects.
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];;
    gkyl_array_release(arrc1->arr);
    gkyl_array_release(arrc2->arr);
  }
  gkyl_free(acs1);
  gkyl_free(acs2);
}

void test_container_pack_accumulate_ho()
{
  // Test an approach to creating an array of array of arrays using the gkyl_container_pack struct.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_containers = 2; // Number of array containers (i.e. no. of arrays).
  int num_packs = 3; // Number of container packs (e.g. number of arrays of containers)

  // Allocate objects.
  struct gkyl_container_pack *cp1 = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  struct gkyl_container_pack *cp2 = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));

  for (int j=0; j<num_packs; j++) {
    cp1[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
    cp2[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));

    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
      arrc1->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
      arrc2->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    }
  }

  // Assign arrays.
  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
      double *arr1_d = arrc1->arr->data;
      for (unsigned i=0; i<arrc1->arr->size; ++i) {
        arr1_d[i] = k*100.0 + i*1.0;
      }

      double *arr2_d  = arrc2->arr->data;
      for (unsigned i=0; i<arrc2->arr->size; ++i) {
        arr2_d[i] = k*200.0 + i*2.0;
      }
    }
  }

  // Accumulate arrays.
  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
      gkyl_array_accumulate(arrc1->arr, 0.5, arrc2->arr);
    }
  }

  // Check results.
  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k];
      double *arr1_d  = arrc1->arr->data;
      for (unsigned i=0; i<arrc1->arr->size; ++i)
        TEST_CHECK( gkyl_compare(arr1_d[i], 2.0*(k*100.0+i*1.0), 1e-14) );
    }
  }

  // Free objects.
  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];;
      gkyl_array_release(arrc1->arr);
      gkyl_array_release(arrc2->arr);
    }
    gkyl_free(acs1);
    gkyl_free(acs2);
  }
  gkyl_free(cp1);
  gkyl_free(cp2);
}

void test_array_bag_accumulate_ho()
{
  // Test an approach to creating an array of array of arrays using the gkyl_array_bag struct.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_arrays = 2; // Number of arrays.
  int num_bags = 3; // Number of array bags (e.g. number of arrays of arrays)

  // Allocate objects.
  struct gkyl_array_bag *ab1 = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  struct gkyl_array_bag *ab2 = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));

  for (int j=0; j<num_bags; j++) {
    ab1[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));
    ab2[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));

    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      innerbag1->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
      innerbag2->arr = gkyl_array_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    }
  }

  // Assign arrays.
  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      struct gkyl_array *arr1 = innerbag1->arr, *arr2 = innerbag2->arr;
      double *arr1_d = arr1->data;
      for (unsigned i=0; i<arr1->size; ++i) {
        arr1_d[i] = k*100.0 + i*1.0;
      }

      double *arr2_d  = arr2->data;
      for (unsigned i=0; i<arr2->size; ++i) {
        arr2_d[i] = k*200.0 + i*2.0;
      }
    }
  }

  // Accumulate arrays.
  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      struct gkyl_array *arr1 = innerbag1->arr, *arr2 = innerbag2->arr;
      gkyl_array_accumulate(arr1, 0.5, arr2);
    }
  }

  // Check results.
  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      struct gkyl_array *arr1 = innerbag1->arr, *arr2 = innerbag2->arr;
      double *arr1_d = arr1->data;
      for (unsigned i=0; i<arr1->size; ++i)
        TEST_CHECK( gkyl_compare(arr1_d[i], 2.0*(k*100.0+i*1.0), 1e-14) );
    }
  }

  // Free objects.
  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      gkyl_array_release(innerbag1->arr);
      gkyl_array_release(innerbag2->arr);
    }
    gkyl_free(ab1[j].bag);
    gkyl_free(ab2[j].bag);
  }
  gkyl_free(ab1);
  gkyl_free(ab2);
}

#ifdef GKYL_HAVE_CUDA

/* Function signatures of kernel calls */
void test_array_container_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1, struct gkyl_array_container *acs2);

void test_array_container_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1, double a, struct gkyl_array_container *acs2);

int test_array_container_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_containers,
  struct gkyl_array_container *acs1);

void test_array_bag_accumulate_dev_assign_cu(int arr_ncomp, int arr_size, int num_bags,
  struct gkyl_array_bag *bag1, struct gkyl_array_bag *bag2);

void test_array_bag_accumulate_dev_accumulate_cu(int arr_ncomp, int arr_size, int num_bags,
  struct gkyl_array_bag *bag1, double a, struct gkyl_array_bag *bag2);

int test_array_bag_accumulate_dev_check_cu(int arr_ncomp, int arr_size, int num_bags,
  struct gkyl_array_bag *bag1);

void test_array_container_accumulate_dev()
{
  // Test an approach to creating an array of arrays using the gkyl_array_container struct.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_containers = 2; // Number of array containers (i.e. no. of arrays).

  // Allocate objects.
  // These hold the host memory pointers.
  struct gkyl_array_container *acs1_ho = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
  struct gkyl_array_container *acs2_ho = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1_ho[k], *arrc2 = &acs2_ho[k];
    arrc1->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    arrc2->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
  }

  // These hold the device-memory pointers.
  struct gkyl_array_container *acs1_dev = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
  struct gkyl_array_container *acs2_dev = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1_ho = &acs1_ho[k], *arrc2_ho = &acs2_ho[k];
    struct gkyl_array_container *arrc1_dev = &acs1_dev[k], *arrc2_dev = &acs2_dev[k];
    arrc1_dev->arr = arrc1_ho->arr->on_dev;
    arrc2_dev->arr = arrc2_ho->arr->on_dev;
  }

  // These are pointers to device memory, and they hold the device-memory pointers.
  struct gkyl_array_container *acs1 = gkyl_cu_malloc(num_containers * sizeof(struct gkyl_array_container));
  struct gkyl_array_container *acs2 = gkyl_cu_malloc(num_containers * sizeof(struct gkyl_array_container));
  gkyl_cu_memcpy(acs1, acs1_dev, num_containers * sizeof(struct gkyl_array_container), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(acs2, acs2_dev, num_containers * sizeof(struct gkyl_array_container), GKYL_CU_MEMCPY_H2D);
  // We can free the _dev ones because we don't need them anymore.
  gkyl_free(acs1_dev);
  gkyl_free(acs2_dev);

  // Assign arrays.
  test_array_container_accumulate_dev_assign_cu(arr_ncomp, arr_size, num_containers, acs1, acs2);

  // Accumulate arrays.
  test_array_container_accumulate_dev_accumulate_cu(arr_ncomp, arr_size, num_containers, acs1, 0.5, acs2);

  // Check results.
  int nfail = test_array_container_accumulate_dev_check_cu(arr_ncomp, arr_size, num_containers, acs1);
  TEST_CHECK( nfail == 0 );

  // Free objects.
  // Note that when you call array_release here you also free the
  // memory that the pointers in acs1/acs2 point to.
  gkyl_cu_free(acs1);
  gkyl_cu_free(acs2);
  for (int k=0; k<num_containers; k++) {
    struct gkyl_array_container *arrc1 = &acs1_ho[k], *arrc2 = &acs2_ho[k];;
    gkyl_array_release(arrc1->arr);
    gkyl_array_release(arrc2->arr);
  }
  gkyl_free(acs1_ho);
  gkyl_free(acs2_ho);
}

void test_container_pack_accumulate_dev()
{
  // Test an approach to creating an array of array of arrays using the gkyl_container_pack struct.
  // MF 2024/10/21: In this test we assume that the container_pack holds pointers to host
  // memory. I think it could be made so that it points to device memory instead too.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_containers = 2; // Number of array containers (i.e. no. of arrays).
  int num_packs = 3; // Number of container packs (e.g. number of arrays of containers)

  // Allocate objects.
  // These hold the host memory pointers.
  struct gkyl_container_pack *cp1_ho = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  struct gkyl_container_pack *cp2_ho = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));

  for (int j=0; j<num_packs; j++) {
    cp1_ho[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
    cp2_ho[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));

    struct gkyl_array_container *acs1 = cp1_ho[j].ac, *acs2 = cp2_ho[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];
      arrc1->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
      arrc2->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    }
  }

  // These hold the device-memory pointers.
  struct gkyl_container_pack *cp1_dev = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  struct gkyl_container_pack *cp2_dev = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  for (int j=0; j<num_packs; j++) {
    cp1_dev[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));
    cp2_dev[j].ac = gkyl_malloc(num_containers * sizeof(struct gkyl_array_container));

    struct gkyl_array_container *acs1_ho = cp1_ho[j].ac, *acs2_ho = cp2_ho[j].ac;
    struct gkyl_array_container *acs1_dev = cp1_dev[j].ac, *acs2_dev = cp2_dev[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1_ho = &acs1_ho[k], *arrc2_ho = &acs2_ho[k];
      struct gkyl_array_container *arrc1_dev = &acs1_dev[k], *arrc2_dev = &acs2_dev[k];
      arrc1_dev->arr = arrc1_ho->arr->on_dev;
      arrc2_dev->arr = arrc2_ho->arr->on_dev;
    }
  }

  // These are pointers to host memory, and they hold the device-memory pointers.
  struct gkyl_container_pack *cp1 = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  struct gkyl_container_pack *cp2 = gkyl_malloc(num_packs * sizeof(struct gkyl_container_pack));
  for (int j=0; j<num_packs; j++) {
    cp1[j].ac = gkyl_cu_malloc(num_containers * sizeof(struct gkyl_array_container));
    cp2[j].ac = gkyl_cu_malloc(num_containers * sizeof(struct gkyl_array_container));
    gkyl_cu_memcpy(cp1[j].ac, cp1_dev[j].ac, num_containers * sizeof(struct gkyl_array_container), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(cp2[j].ac, cp2_dev[j].ac, num_containers * sizeof(struct gkyl_array_container), GKYL_CU_MEMCPY_H2D);
  }
  // We can free the _dev ones because we don't need them anymore.
  for (int j=0; j<num_packs; j++) {
    gkyl_free(cp1_dev[j].ac);
    gkyl_free(cp2_dev[j].ac);
  }
  gkyl_free(cp1_dev);
  gkyl_free(cp2_dev);

  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1[j].ac, *acs2 = cp2[j].ac;

    // Assign arrays.
    test_array_container_accumulate_dev_assign_cu(arr_ncomp, arr_size, num_containers, acs1, acs2);

    // Accumulate arrays.
    test_array_container_accumulate_dev_accumulate_cu(arr_ncomp, arr_size, num_containers, acs1, 0.5, acs2);

    // Check results.
    int nfail = test_array_container_accumulate_dev_check_cu(arr_ncomp, arr_size, num_containers, acs1);
    TEST_CHECK( nfail == 0 );
  }

  // Free objects.
  // Note that when you call array_release here you also free the
  // memory that the pointers in cp1/cp2 point to.
  for (int j=0; j<num_packs; j++) {
    gkyl_cu_free(cp1[j].ac);
    gkyl_cu_free(cp2[j].ac);
  }
  gkyl_free(cp1);
  gkyl_free(cp2);
  for (int j=0; j<num_packs; j++) {
    struct gkyl_array_container *acs1 = cp1_ho[j].ac, *acs2 = cp2_ho[j].ac;
    for (int k=0; k<num_containers; k++) {
      struct gkyl_array_container *arrc1 = &acs1[k], *arrc2 = &acs2[k];;
      gkyl_array_release(arrc1->arr);
      gkyl_array_release(arrc2->arr);
    }
    gkyl_free(acs1);
    gkyl_free(acs2);
  }
  gkyl_free(cp1_ho);
  gkyl_free(cp2_ho);
}

void test_array_bag_accumulate_dev()
{
  // Test an approach to creating an array of array of arrays using the gkyl_array_bag struct.
  // MF 2024/10/21: In this test we assume that (outer) the array_bag holds pointers to host
  // memory. I think it could be made so that it points to device memory instead too.
  int arr_ncomp = 1; // Number of components of each array.
  int arr_size = 10; // Number of elements/cells of each array.
  int num_arrays = 2; // Number of array containers (i.e. no. of arrays).
  int num_bags = 3; // Number of container packs (e.g. number of arrays of containers)

  // Allocate objects.
  // These hold the host memory pointers.
  struct gkyl_array_bag *ab1_ho = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  struct gkyl_array_bag *ab2_ho = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  for (int j=0; j<num_bags; j++) {
    ab1_ho[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));
    ab2_ho[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));

    struct gkyl_array_bag *bag1 = &ab1_ho[j], *bag2 = &ab2_ho[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      innerbag1->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
      innerbag2->arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ncomp, arr_size);
    }
  }

  // These hold the device-memory pointers.
  struct gkyl_array_bag *ab1_dev = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  struct gkyl_array_bag *ab2_dev = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  for (int j=0; j<num_bags; j++) {
    ab1_dev[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));
    ab2_dev[j].bag = gkyl_malloc(num_arrays * sizeof(struct gkyl_array_bag));

    struct gkyl_array_bag *bag1_ho = &ab1_ho[j], *bag2_ho = &ab2_ho[j];
    struct gkyl_array_bag *bag1_dev = &ab1_dev[j], *bag2_dev = &ab2_dev[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1_ho = &bag1_ho->bag[k], *innerbag2_ho = &bag2_ho->bag[k];
      struct gkyl_array_bag *innerbag1_dev = &bag1_dev->bag[k], *innerbag2_dev = &bag2_dev->bag[k];
      innerbag1_dev->arr = innerbag1_ho->arr->on_dev;
      innerbag2_dev->arr = innerbag2_ho->arr->on_dev;
    }
  }

  // These are pointers to host memory, and they hold the device-memory pointers.
  struct gkyl_array_bag *ab1 = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  struct gkyl_array_bag *ab2 = gkyl_malloc(num_bags * sizeof(struct gkyl_array_bag));
  for (int j=0; j<num_bags; j++) {
    ab1[j].bag = gkyl_cu_malloc(num_arrays * sizeof(struct gkyl_array_bag));
    ab2[j].bag = gkyl_cu_malloc(num_arrays * sizeof(struct gkyl_array_bag));

    gkyl_cu_memcpy(ab1[j].bag, ab1_dev[j].bag, num_arrays * sizeof(struct gkyl_array_bag), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(ab2[j].bag, ab2_dev[j].bag, num_arrays * sizeof(struct gkyl_array_bag), GKYL_CU_MEMCPY_H2D);
  }
  // We can free the _dev ones because we don't need them anymore.
  for (int j=0; j<num_bags; j++) {
    gkyl_free(ab1_dev[j].bag);
    gkyl_free(ab2_dev[j].bag);
  }
  gkyl_free(ab1_dev);
  gkyl_free(ab2_dev);

  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1[j], *bag2 = &ab2[j];

    // Assign arrays.
    test_array_bag_accumulate_dev_assign_cu(arr_ncomp, arr_size, num_arrays, bag1->bag, bag2->bag);

    // Accumulate arrays.
    test_array_bag_accumulate_dev_accumulate_cu(arr_ncomp, arr_size, num_arrays, bag1->bag, 0.5, bag2->bag);

    // Check results.
    int nfail = test_array_bag_accumulate_dev_check_cu(arr_ncomp, arr_size, num_arrays, bag1->bag);
    TEST_CHECK( nfail == 0 );
  }

  // Free objects.
  // Note that when you call array_release here you also free the
  // memory that the pointers in ab1/ab2 point to.
  for (int j=0; j<num_bags; j++) {
    gkyl_cu_free(ab1[j].bag);
    gkyl_cu_free(ab2[j].bag);
  }
  gkyl_free(ab1);
  gkyl_free(ab2);
  for (int j=0; j<num_bags; j++) {
    struct gkyl_array_bag *bag1 = &ab1_ho[j], *bag2 = &ab2_ho[j];
    for (int k=0; k<num_arrays; k++) {
      struct gkyl_array_bag *innerbag1 = &bag1->bag[k], *innerbag2 = &bag2->bag[k];
      gkyl_array_release(innerbag1->arr);
      gkyl_array_release(innerbag2->arr);
    }
    gkyl_free(ab1_ho[j].bag);
    gkyl_free(ab2_ho[j].bag);
  }
  gkyl_free(ab1_ho);
  gkyl_free(ab2_ho);
}

#endif

TEST_LIST = {
  { "array_container_accumulate_ho", test_array_container_accumulate_ho },
  { "container_pack_accumulate_ho", test_container_pack_accumulate_ho },
  { "array_bag_accumulate_ho", test_array_bag_accumulate_ho },
#ifdef GKYL_HAVE_CUDA
  { "array_container_accumulate_dev", test_array_container_accumulate_dev },
  { "container_pack_accumulate_dev", test_container_pack_accumulate_dev },
  { "array_bag_accumulate_dev", test_array_bag_accumulate_dev },
#endif
  { NULL, NULL },
};

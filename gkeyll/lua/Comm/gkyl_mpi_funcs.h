// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <gkyl_mpi_macros.h>

// Sizes of various objects
DECL_GET_MPI_OBJ_SIZE(MPI_Status);
DECL_GET_MPI_OBJ_SIZE(MPI_Aint);
DECL_GET_MPI_OBJ_SIZE(MPI_Request);
DECL_GET_MPI_OBJ_PTR_SIZE(MPI_Status);
DECL_GET_MPI_OBJ_PTR_SIZE(MPI_Request);
  
// Pre-defined objects and constants
DECL_GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);
DECL_GET_MPI_OBJECT(Comm, MPI_COMM_NULL);
DECL_GET_MPI_OBJECT(Comm, MPI_COMM_SELF);
DECL_GET_MPI_OBJECT(Info, MPI_INFO_NULL);
  
DECL_GET_MPI_OBJECT(Request, MPI_REQUEST_NULL);
DECL_GET_MPI_OBJECT_PTR(Status, MPI_STATUS_IGNORE);
DECL_INT_OBJECT(MPI_PROC_NULL);

// Datatypes
DECL_GET_MPI_OBJECT(Datatype, MPI_C_BOOL);
DECL_GET_MPI_OBJECT(Datatype, MPI_CHAR);
DECL_GET_MPI_OBJECT(Datatype, MPI_BYTE);
DECL_GET_MPI_OBJECT(Datatype, MPI_SHORT);
DECL_GET_MPI_OBJECT(Datatype, MPI_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_LONG);
DECL_GET_MPI_OBJECT(Datatype, MPI_FLOAT);
DECL_GET_MPI_OBJECT(Datatype, MPI_DOUBLE);
DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_CHAR);
DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_SHORT);
DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED);
DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_LONG);
DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE);
DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_LONG_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_FLOAT_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_DOUBLE_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_SHORT_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_2INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE_INT);
DECL_GET_MPI_OBJECT(Datatype, MPI_PACKED);

// Ops
DECL_GET_MPI_OBJECT(Op, MPI_MAX);
DECL_GET_MPI_OBJECT(Op, MPI_MIN);
DECL_GET_MPI_OBJECT(Op, MPI_SUM);
DECL_GET_MPI_OBJECT(Op, MPI_PROD);
DECL_GET_MPI_OBJECT(Op, MPI_LAND);
DECL_GET_MPI_OBJECT(Op, MPI_BAND);
DECL_GET_MPI_OBJECT(Op, MPI_LOR);
DECL_GET_MPI_OBJECT(Op, MPI_BOR);
DECL_GET_MPI_OBJECT(Op, MPI_LXOR);
DECL_GET_MPI_OBJECT(Op, MPI_BXOR);
DECL_GET_MPI_OBJECT(Op, MPI_MINLOC);
DECL_GET_MPI_OBJECT(Op, MPI_MAXLOC);

// error codes
DECL_INT_OBJECT(MPI_SUCCESS);

// Constants
DECL_INT_OBJECT(MPI_COMM_TYPE_SHARED);
DECL_INT_OBJECT(MPI_UNDEFINED);
DECL_INT_OBJECT(MPI_ORDER_C);
DECL_INT_OBJECT(MPI_ORDER_FORTRAN);
DECL_INT_OBJECT_PTR(MPI_UNWEIGHTED);
DECL_VOID_OBJECT_PTR(MPI_BOTTOM);

// Some utility functions to allow accessing non-opaque MPI types

// Gkyl structs holding status and requests.
typedef struct {
  MPI_Request *req;
} gkyl_MPI_Request;

typedef struct {
  MPI_Status *stat;
} gkyl_MPI_Status;

typedef struct {
  MPI_Request *req;
  MPI_Status *stat;
} gkyl_MPI_Request_Status;

// Functions to allocate and free structs holding requests and statuses.
void gkyl_MPI_Request_alloc(gkyl_MPI_Request *rs, int num);
void gkyl_MPI_Request_release(gkyl_MPI_Request *rs);
void gkyl_MPI_Status_alloc(gkyl_MPI_Status *ss, int num);
void gkyl_MPI_Status_release(gkyl_MPI_Status *ss);
void gkyl_MPI_Request_Status_alloc(gkyl_MPI_Request_Status *rss, int num);
void gkyl_MPI_Request_Status_release(gkyl_MPI_Request_Status *rss);

// Functions to fetch members of status.
int gkyl_mpi_get_status_SOURCE(const MPI_Status* instat, int off);
int gkyl_mpi_get_status_TAG(const MPI_Status* instat, int off);
int gkyl_mpi_get_status_ERROR(const MPI_Status* instat, int off);

// Get count from a status (which may be one of several in an array of
// statuses).
int gkyl_mpi_get_status_count(const MPI_Status *instat, MPI_Datatype datatype, int *count, int off);

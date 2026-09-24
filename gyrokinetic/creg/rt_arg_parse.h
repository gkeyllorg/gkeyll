#pragma once

#include <gkyl_basis.h>
#include <gkyl_mp_scheme.h>

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define APP_ARGS_DEFAULT_FILE_NAME "name-not-set"

#define APP_ARGS_CHOOSE(val, def) ((val) > 0 ? (val) : (def))

struct gkyl_app_args {
  bool use_gpu; // should this be run on GPU?
  bool use_mpi; // should this be run on MPI?  
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  bool trace_mem; // should we trace memory allocation/deallocations?
  int num_steps; // number of steps
  int num_threads; // number of threads
  int xcells[3]; // configuration space cells
  int vcells[3]; // velocity space cells
  int cuts[3]; // domain decomposition "cuts"
  char file_name[1024]; // name of input file
  char app_name[128]; // basename of argv[0] — use with snprintf/strcpy to set char[] name fields (cannot assign directly in compound literals)
  enum gkyl_basis_type basis_type; // type of basis functions to use
  enum gkyl_mp_recon mp_recon; // the XX in MP-XX
  bool skip_limiters; // should we skip limiters?
  bool is_restart; // is this a restarted sim?
  int restart_frame; // frame to restart from
  char opt_args[128]; // optional arguments
};

static int
get_basis_type(const char *nm)
{
  if (strcmp(nm, "ms") == 0) {
    return GKYL_BASIS_MODAL_SERENDIPITY;
  }
  else if (strcmp(nm, "mt") == 0) {
    return GKYL_BASIS_MODAL_TENSOR;
  }
  return -1;
}

static int
get_mp_recon_type(const char *nm)
{
  if (strcmp(nm, "u1") == 0) {
    return GKYL_MP_U1;
  }
  else if (strcmp(nm, "u3") == 0) {
    return GKYL_MP_U3;
  }
  else if (strcmp(nm, "u5") == 0) {
    return GKYL_MP_U5;
  }
  else if (strcmp(nm, "c2") == 0) {
    return GKYL_MP_C2;
  }
  else if (strcmp(nm, "c4") == 0) {
    return GKYL_MP_C4;
  }
  else if (strcmp(nm, "c6") == 0) {
    return GKYL_MP_C6;
  }  
  
  return -1;
}

static struct gkyl_app_args
parse_app_args(int argc, char **argv)
{
  bool use_gpu = false;
  bool use_mpi = false;
  bool step_mode = false;
  bool trace_mem = false;
  bool is_restart = false;
  int restart_frame = INT_MAX;
  bool skip_limiters = false;
  int num_steps = INT_MAX;
  int num_threads = 1; // by default use only 1 thread

  struct gkyl_app_args args = {
    .xcells = { 0 },
    .vcells = { 0 },
    .cuts = { 1, 1, 1 },
  };

  strcpy(args.file_name, APP_ARGS_DEFAULT_FILE_NAME); // default
  args.basis_type = GKYL_BASIS_MODAL_SERENDIPITY;

  int c;
  while ((c = getopt(argc, argv, "+hjgmMt:s:i:b:x:y:z:u:v:w:r:c:d:e:o:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: <app_name> -g -m -s nsteps -t nthreads -i inp -b [ms|mt] -x NX -y NY -z NZ -u VX -v VY -w VZ\n");
        printf(" All flags and parameters are optional.\n");
        printf(" -g     Run on GPUs if GPUs are present and code built for GPUs\n");
        printf(" -M     Run with MPI if code built with MPI\n");        
        printf(" -sN    Only run N steps of simulation\n");
        printf(" -tN    Use N threads (when available)\n");
        printf(" -b     Basis function to use (ms: Modal serendipity; mt: Modal tensor-product)\n");
        printf("        (Ignored for finite-volume solvers)\n");
        printf(" -j     Recovery scheme. One of u1, u3, u5, c2, c4, c6\n");
        printf("        (Only used for MP-XX solvers)\n");
        printf(" -l     Turn off limiters\n");
        printf(" -rN    Restart the simulation from frame N\n");
        printf(" -m     Turn on memory allocation/deallocation tracing\n");
        printf(" -o     Optional arguments (as string, requires parsing)\n");
        printf("\n");
        printf(" Grid resolution in configuration space:\n");
        printf(" -xNX -yNY -zNZ\n");
        printf(" Grid resolution in velocity space:\n");
        printf(" -uVX -vVY -wVZ\n");
        printf(" Domain decomposition in each direction:\n");
        printf(" -cPX -dPY -ePZ\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;

      case 'M':
        use_mpi = true;
        break;        

      case 'm':
        trace_mem = true;
        break;

      case 'l':
        skip_limiters = true;
        break;

      case 'r':
        is_restart = true;
        restart_frame = atoi(optarg);
        break;        
      
      case 's':
        step_mode = true;
        num_steps = atoi(optarg);
        break;

      case 't':
        num_threads = atoi(optarg);
        break;

      case 'c':
        args.cuts[0] = atoi(optarg);
        break;
        
      case 'd':
        args.cuts[1] = atoi(optarg);
        break;
        
      case 'e':
        args.cuts[2] = atoi(optarg);
        break;        

      case 'x':
        args.xcells[0] = atoi(optarg);
        break;
        
      case 'y':
        args.xcells[1] = atoi(optarg);
        break;
        
      case 'z':
        args.xcells[2] = atoi(optarg);
        break;

      case 'u':
        args.vcells[0] = atoi(optarg);
        break;
        
      case 'v':
        args.vcells[1] = atoi(optarg);
        break;
        
      case 'w':
        args.vcells[2] = atoi(optarg);
        break;        

      case 'i':
        strcpy(args.file_name, optarg);
        break;

      case 'b':
        args.basis_type = get_basis_type(optarg);
        assert(args.basis_type != -1);
        break;

     case 'j':
        args.mp_recon = get_mp_recon_type(optarg);
        assert(args.mp_recon != -1);
        break;        

      case 'o':
        assert(strlen(optarg) < sizeof(args.opt_args));
        strcpy(args.opt_args, optarg);
        break;

      case '?':
        break;
    }
  }
  
  args.use_gpu = use_gpu;
  args.use_mpi = use_mpi;
  args.trace_mem = trace_mem;
  args.step_mode = step_mode;
  args.num_steps = num_steps;
  args.num_threads = num_threads;
  args.skip_limiters = skip_limiters;
  args.is_restart = is_restart;
  args.restart_frame = restart_frame;

  // Derive app_name from the executable name (basename of argv[0]).
  if (argc > 0 && argv[0]) {
    const char *slash = strrchr(argv[0], '/');
    const char *base = slash ? slash + 1 : argv[0];
    strncpy(args.app_name, base, sizeof(args.app_name) - 1);
    args.app_name[sizeof(args.app_name) - 1] = '\0';
  } else {
    strncpy(args.app_name, "unknown", sizeof(args.app_name) - 1);
  }

  return args;
}

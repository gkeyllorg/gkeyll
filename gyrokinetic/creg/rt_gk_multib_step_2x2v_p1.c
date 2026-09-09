#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

void shaped_pfunc_lower_outer(double s, double* RZ){
  double p0[2] = {5.488-0.6,-8.600};
  double p1[2] = {5.855-0.6,-8.52318};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  double p0[2] = {5.488-0.6,8.600};
  double p1[2] = {5.855-0.6,8.52318};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}

//old inner plates
void shaped_pfunc_upper_inner(double s, double* RZ){
  RZ[0] = 1.651 + (1.8 - 1.651)*s;
  RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
  RZ[0] = 1.651 + (1.8 - 1.651)*s;
  RZ[1] = -(6.331 + (6.777 - 6.331)*s);
}

struct gkyl_gk_block_geom*
create_gk_block_geom(void)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 12);

  /* Block layout and coordinates

   x  
   ^  
   |
   4  +------------------+------------------+------------------+
   |  |b1                |b2                |b3                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   3  +------------------+------------------+------------------+
   |  |b0               x|o b10            %|$ b4              |
   |  |lower outer PF   x|o outer core     %|$ upper outer PF  |
   |  |                 x|o                %|$                 |
   |  +------------------+------------------+------------------+
   2  +------------------+------------------+------------------+
   |  |b9               x|o b11            %|$ b5              |
   |  |lower inner PF   x|o inner core     %|$ upper inner PF  |
   |  |                 x|o                %|$                 |
   1  +------------------+------------------+------------------+
   |  |b8                |b7                |b6                |
   |  |lower inner SOL   |middle inner SOL  |upper inner SOL   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+

      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);

  double wout = 0.1;
  double win = 0.05;
  double wcore = 0.1;
  double wpf = 0.05;

  double psi_lo_outer_sol = psisep - wout;
  double psi_up_outer_sol = psisep;

  double psi_lo_core = psisep;
  double psi_up_core = psisep + wcore;

  double psi_lo_pf = psisep;
  double psi_up_pf = psisep + wpf;

  double psi_lo_inner_sol = psisep - win ;
  double psi_up_inner_sol = psisep;

  int npsi_outer_sol = 4;
  int npsi_core = 4;
  int npsi_pf = 2;
  int npsi_inner_sol = 2;

  double ntheta_lower_inner  = 4;
  double ntheta_middle_inner = 4*2;
  double ntheta_upper_inner  = 4;

  double ntheta_lower_outer = 6;
  double ntheta_middle_outer = 6*2;
  double ntheta_upper_outer = 6;

  double ncuts_lower_inner  = 1;//ntheta_lower_inner;
  double ncuts_middle_inner = 1;//ntheta_middle_inner;
  double ncuts_upper_inner  = 1;//ntheta_upper_inner;

  double ncuts_lower_outer  = 1;//ntheta_lower_outer;
  double ncuts_middle_outer = 1;//ntheta_middle_outer;
  double ncuts_upper_outer  = 1;//ntheta_upper_outer;

  double zinner = 6.34;
  double zouter = 8.29;
  double rright_out = 5.2;

  // Note that for tokamak multi-block simulations, 
  // these theta limits are just placeholders and will be 
  // reset in the multi-block app.
  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

   // block 0. Lower outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_outer },
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_R,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 2. Middle outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psi_up_outer_sol, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle_outer},
      .cuts = { 1, ncuts_middle_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 10, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 3. Upper outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psi_up_outer_sol, theta_up },
      .cells = { npsi_outer_sol, ntheta_upper_outer},
      .cuts = { 1, ncuts_upper_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

  // block 4. Upper outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_outer },
      .cuts = { 1, ncuts_upper_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_UP_R,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // block 5. Upper inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_inner },
      .cuts = { 1, ncuts_upper_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_UP_L,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 6. Upper inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_upper_inner},
      .cuts = { 1, ncuts_upper_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP,
          .rleft = 2.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
   
  // block 7. Middle inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_middle_inner},
      .cuts = { 1, ncuts_middle_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID,
          .rleft = 2.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 11, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 8. Lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 8, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_lower_inner},
      .cuts = { 1, ncuts_lower_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO,
          .rleft = 2.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 9, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_inner },
      .cuts = { 1, ncuts_lower_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_L,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 8, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );


  // block 10. outer core.
  gkyl_gk_block_geom_set_block(bgeom, 10, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_core, theta_lo},
      .upper = { psi_up_core,  theta_up},
      .cells = { npsi_core, ntheta_middle_outer},
      .cuts = { 1, ncuts_middle_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_CORE_R,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 2.0,        // closest R to inboard SOL
          .rmin = 1.58,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 11. Inner Core.
  gkyl_gk_block_geom_set_block(bgeom, 11, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_core, theta_lo },
      .upper = { psi_up_core, theta_up },
      .cells = { npsi_core, ntheta_middle_inner},
      .cuts = { 1, ncuts_middle_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
//        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_CORE_L,
          .rclose = 0.0,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 2.0,        // closest R to inboard SOL
          .rmin = 1.58,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
}
struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double vtIon;
  double vtElc;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double temp_fac;
  double n_fac;
  double cx;
  double cz;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge
  double n_fac = 1.0;
  double temp_fac = 9.0/n_fac;

  double Te = 1500.0/temp_fac*eV;
  double Ti = 1500.0/temp_fac*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);

  // Source parameters
  double nsource = 1.675e22*2.5*temp_fac/6.0;
  double T_source = 3000.0*eV/temp_fac;
  double cx = 0.0065612;
  double cz = 0.4916200*1.4;
  double psisep = 1.5093065418975686; //This is the value from efit

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 8.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 12;

  double t_end = 100.0e-10; 
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .nsource = nsource,
    .T_source = T_source, 
    .temp_fac = temp_fac,
    .n_fac = n_fac,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void
init_density_core(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //1.5e20 at inner core and 2e19 at sep
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 3.0e20;
  double intercept = -4.22792e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_outer(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 2.7e20;
  double intercept = -3.77513e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_inner(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 5.4e20;
  double intercept = -7.85026e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_pf(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = -5.4e20;
  double intercept = 8.45026e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double nsource = app->nsource;
  if(x >= 1.5593065418975687)
    fout[0] = nsource;
  else 
    fout[0] = nsource*1.0e-5;
}

void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Te;
  fout[0] = T;
}

void
init_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Ti;
  fout[0] = T;
}
void
source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->T_source;
  fout[0] = T;
}

static inline void
mapc2p_vel_elc(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_elc = app->mu_max_elc;
  double vpar_max_elc = app->vpar_max_elc;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_elc*2.0*pow(cvpar,2);

  mu = mu_max_elc * (cmu * cmu);

  // Set rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
}

static inline void
mapc2p_vel_ion(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_ion = app->mu_max_ion;
  double vpar_max_ion = app->vpar_max_ion;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_ion*2.0*pow(cvpar,2);

  mu = mu_max_ion * (cmu * cmu);

  // Set rescaled ion velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar ; vp[1] = mu;
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_step_ctx *app = ctx;

  fout[0] = 0.22; // Diffusivity [m^2/s].
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  // construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[12];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[8] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 8,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[9] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 9,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[10] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 10,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

  elc_blocks[11] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 11,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };


  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    //block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    //block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},

  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
    }, 

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },

    .duplicate_across_blocks = false,
    .blocks = elc_blocks,
    .num_physical_bcs = 20,
    .bcs = elc_phys_bcs,
  };

  // Ion Species
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[12];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[8] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 8,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[9] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 9,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

 ion_blocks[10] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 10,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 4,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

 ion_blocks[11] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 11,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 4,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    //block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    //block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},

  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
    }, 
  
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },

    .duplicate_across_blocks = false,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    // block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    // block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 12,
    .bcs = field_phys_bcs,
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_multib *app_inp = gkyl_calloc(1,sizeof(struct gkyl_gyrokinetic_multib));
  app_inp->cdim = ctx.cdim;
  app_inp->poly_order = 1;
  app_inp->basis_type = app_args.basis_type;
  app_inp->use_gpu = app_args.use_gpu;
  app_inp->cfl_frac = 1.0;
  app_inp->gk_block_geom = bgeom;
  app_inp->num_species = 2;
  app_inp->species[0] = elc;
  app_inp->species[1] = ion;
  app_inp->field = field;
  app_inp->comm = comm;

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp->name, sizeof(app_inp->name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = *app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);
  gkyl_free(app_inp);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}

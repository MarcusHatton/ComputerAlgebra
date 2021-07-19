// void regrid_AMPL( const int interp_stencil_size, const REAL AMPL_NEW, const int Nxx_plus_2NGHOSTS[3], REAL *xx[3], REAL *helper_gfs, REAL *in_and_out_gfs ) {

//   /* Useful auxiliary variables */
//   const REAL inv_SINHW          = 1.0/SINHW;
//   const REAL sinh_inv_SINHW     = sinh( inv_SINHW );
//   const REAL inv_sinh_inv_SINHW = 1.0 / sinh_inv_SINHW;

//   for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {

//     /* Find the new value of r after the regrid */
//     const REAL r_star   = AMPL_NEW * sinh( xx[0][i0] * inv_SINHW ) * inv_sinh_inv_SINHW;

//     /* Map the new value of r onto the old x grid */
//     const REAL x_star   = SINHW * asinh( r_star * sinh_inv_SINHW / AMPL );

//     /* Find the index in the x[0] array such that | x[0][idx] - x_star | is minimal */
//     const int idx = bisection_idx_finder(x_star, Nxx_plus_2NGHOSTS[0], xx[0]);

//     /* Set up minimum and maximum interpolation indices */
//     const int idxmin = MAX(0,idx-interp_stencil_size/2-1);
//     const int idxmax = idxmin + interp_stencil_size;

//     /* printf("x_star: %e | x[0][%d]: %e\n",x_star,idxmin,x[0][idxmin]); */

//     /* Compute l_i(x) for the Lagrange polynomial interpolation */
//     REAL l_i_of_x[interp_stencil_size];
//     for(int i=idxmin;i<idxmax;i++) {
//       REAL numer = 1.0;
//       REAL denom = 1.0;
//       for(int j=idxmin;j<i;j++) {
// 	numer *= x_star   - xx[0][j];
// 	denom *= xx[0][i] - xx[0][j];
//       }
//       for(int j=i+1;j<idxmax;j++) {
// 	numer *= x_star   - xx[0][j];
// 	denom *= xx[0][i] - xx[0][j];
//       }
//       l_i_of_x[i-idxmin] = numer/denom;      
//     }

//     /* Perform the interpolation */
//     for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
//       for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
// 	for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
// 	  const int interp_idx   = IDX4S(which_gf,i0,i1,i2);
// 	  helper_gfs[interp_idx] = 0.0;
// 	  for( int i=idxmin; i<idxmax; i++ ) {
// 	    helper_gfs[interp_idx] += l_i_of_x[i-idxmin] * in_and_out_gfs[IDX4S(which_gf,i,i1,i2)];
// 	  }
// 	}
//       }
//     }

//   } // END OF for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++)

//   /* Now update the gridfunctions */
//   for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
//     for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
//       for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
// 	for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {
// 	  const int index       = IDX4S(which_gf,i0,i1,i2);
// 	  in_and_out_gfs[index] = helper_gfs[index];
// 	}
//       }
//     }
//   }

//   /* Update the amplitude */
//   AMPL = AMPL_NEW;

// }

void regrid_SINHW( const paramstruct *restrict params,
                   const int interp_stencil_size,
                   const REAL SINHW_NEW,
                   REAL *restrict xx[3],
                   REAL *restrict helper_gfs,
                   REAL *restrict in_and_out_gfs ) {
#include "set_Cparameters.h"
  
  /* Useful auxiliary variables */
  const REAL inv_SINHW          = 1.0/SINHW;
  const REAL sinh_inv_SINHW     = sinh( inv_SINHW );
  const REAL inv_sinh_inv_SINHW = 1.0 / sinh_inv_SINHW;

  const REAL inv_SINHW_NEW          = 1.0/SINHW_NEW;
  const REAL sinh_inv_SINHW_NEW     = sinh( inv_SINHW_NEW );
  const REAL inv_sinh_inv_SINHW_NEW = 1.0 / sinh_inv_SINHW_NEW;

  for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS0-NGHOSTS;i0++) {

    /* Find the new value of r after the regrid */
    const REAL r_star   = AMPL * sinh( xx[0][i0] * inv_SINHW_NEW ) * inv_sinh_inv_SINHW_NEW;

    /* Map the new value of r onto the old x grid */
    const REAL x_star   = SINHW * asinh( r_star * sinh_inv_SINHW / AMPL );

    /* Find the index in the x[0] array such that | x[0][idx] - x_star | is minimal */
    const int idx = bisection_idx_finder(x_star, Nxx_plus_2NGHOSTS0, xx[0]);

    /* Set up minimum and maximum interpolation indices */
    const int idxmin = MAX(0,idx-interp_stencil_size/2-1);
    const int idxmax = idxmin + interp_stencil_size;

    /* printf("x_star: %e | x[0][%d]: %e\n",x_star,idxmin,x[0][idxmin]); */

    /* Compute l_i(x) for the Lagrange polynomial interpolation */
    REAL l_i_of_x[interp_stencil_size];
    for(int i=idxmin;i<idxmax;i++) {
      REAL numer = 1.0;
      REAL denom = 1.0;
      for(int j=idxmin;j<i;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      for(int j=i+1;j<idxmax;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      l_i_of_x[i-idxmin] = numer/denom;      
    }

    /* Perform the interpolation */
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
      for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS2-NGHOSTS;i2++) {
	for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS1-NGHOSTS;i1++) {
	  const int interp_idx   = IDX4S(which_gf,i0,i1,i2);
	  helper_gfs[interp_idx] = 0.0;
	  for( int i=idxmin; i<idxmax; i++ ) {
	    helper_gfs[interp_idx] += l_i_of_x[i-idxmin] * in_and_out_gfs[IDX4S(which_gf,i,i1,i2)];
	  }
	}
      }
    }

  } // END OF for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++)

  /* Now update the gridfunctions */
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS2-NGHOSTS;i2++) {
      for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS1-NGHOSTS;i1++) {
	for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS0-NGHOSTS;i0++) {
	  const int index       = IDX4S(which_gf,i0,i1,i2);
	  in_and_out_gfs[index] = helper_gfs[index];
	}
      }
    }
  }

}

void regrid(const int n,
            const REAL t,
            const REAL sinhW_new,
            const bc_struct *restrict bcstruct,
            const rfm_struct *restrict rfmstruct,
            paramstruct *restrict params,
            int *restrict N_final,
            REAL *restrict t_final,
            REAL *restrict dt,
            REAL *restrict xx[3],
            REAL *restrict aux_gfs,
            REAL *restrict in_gfs) {

  // Set interpolation stencil size
  const int regrid_stencil_size = 2*NGHOSTS - 1;

  // Perform the regrid
  regrid_SINHW(params,regrid_stencil_size,sinhW_new,xx,aux_gfs,in_gfs);

  // Update SINHW
  params->SINHW = sinhW_new;

  // Recompute rfm_struct
#include "set_Cparameters.h"
#include "rfm_files/rfm_struct__define-pointer.h"

  // Apply boundary conditions to the evolved gridfunctions
  // and enforce the deteterminant of gammabar constraint
  apply_bcs_curvilinear(params,bcstruct,NUM_EVOL_GFS,evol_gf_parity,in_gfs);
  enforce_detgammabar_constraint(rfmstruct,params,in_gfs);

  // Update dt
  *dt      = find_timestep(params, xx);

  // Update t_final
  *t_final = MIN(*t_final,t+params->AMPL);

  // Update N_final
  int remaining_iterations = (int)( ( (*t_final) - t ) / (*dt) + 0.5);
  *N_final = n + remaining_iterations;

}

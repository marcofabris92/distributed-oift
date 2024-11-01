/*
 * cost.c
 *
 *    incremental cost - L_2 traj error (squared)
 *
 *      l(x,u,wct)   [ w_c(t) = x_des ]
 *            = nAg*|| xB - xB_des ||^2_QB / 2 + || u ||^2_RD / 2 + kF*FF(p) + kA*AA(dp)
 *      m(x,u,wct)
 *						= || xB - xB_des ||^2_P1
 *
 *
 *		FF(x) = 1/2 * sum_i=1^nAg sum_for all j\neq i [sigma(||pi-pj||^2)]
 *
 *
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),  Q(8),Qsafe(64), S(16), R(32)
 *
 * state             (x)
 *
 *   p1, ... , pn, dp1, ... , dpn
 *
 * input             (u)
 *
 *   u1, ..., un
 *
 * exogenous input   (wct)
 *
 *   xB_des (pB, dpB)
 *
 * Marco Fabris, Boulder, 06/29/2018
 * checked
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),  Q(8),Qsafe(64), S(16), R(32)
 */

#include <math.h>
#include "sys_sizes.h"
#include "cost_params.h"
#include "sigma.c"


/*
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),  Q(8),Qsafe(64), S(16), R(32)
*/
void
 cost(
      double *x, double *u, double *wct,
      int ders,
      double *lxu,
      double *lxu_x_, double *lxu_u_,
      double *lxu_x_x_, double *lxu_x_u_, double *lxu_u_u_){

#define NX  (NS)
#define NU  (NI)

#define lxu_x_x(i,j)      lxu_x_x_[ ((i)-1) + ((j)-1)*(NX) ]  	  // columnwise
#define lxu_x_u(i,j)      lxu_x_u_[ ((i)-1) + ((j)-1)*(NX) ]  	  // columnwise
#define lxu_u_u(i,j)      lxu_u_u_[ ((i)-1) + ((j)-1)*(NU) ]  	  // columnwise

#define nAg2  (nAg*nAg)          // number of agents squared
#define INTS  ((nAg*(nAg-1))/2)	 // number of interactions between different agents
#define SQRT2 (sqrt(2.0))
#define SQRT3 (sqrt(3.0))
#define RADIUS5 (1.0/(2.0*sin(M_PI/5.0)))
#define DIAG5   (sqrt(2.0-2.0*cos(3.0*M_PI/5.0)))
#define DIAGP  (2*sin(M_PI*54/180))

  double *xB_des = wct;
  double dx[NX], du[NU], dxB[2*DIM];

  int do_l, do_a, do_b, do_Q, do_S, do_R, do_Qsafe;

  double l_st, l_in, l_fo;

  // determine what jobs to do
  /*
   *    ders - which ders? use binary code to specify
   *    1 - lxu
   *    2 - lxu_x = a
   *    4 - lxu_u = b
   *    8 - lxu_x_x = Q
   *   16 - lxu_x_u = S
   *   32 - lxu_u_u = R
   */
  do_l     =  ders & 1       ;
  do_a     = (ders & 2)  >> 1;
  do_b     = (ders & 4)  >> 2;
  do_Q     = (ders & 8)  >> 3;
  do_S     = (ders & 16) >> 4;
  do_R     = (ders & 32) >> 5;
  do_Qsafe = (ders & 64) >> 6;

  // a = C'*QB*(xB-xB_des) + nabla_x_FF    			| NS
  // b = RD*u   						    								| NI
  // Q = C'*QB*C + H_x_x_FF 		        				| NSxNS
  // R = RD								    									| NIxNI
  // S = 0 								    									| NSxNI


  // C = 1/nAg*[I_D ... I_D 0_D ... 0_D
  //            0_D ... 0_D I_D ... I_D]

  // gradient and Hessian of kF*FF wrt x
  // gradient_x_FF  = kF*[gradient_p_p_FF         =  [gradient_p_p_FF
  //                      gradient_pdot_pdot_FF]           0          ]
  // hessian_x_x_FF = kF*[hessian_p_p_FF       | hessian_p_pdot_FF          = [ hessian_p_p_FF  |  0
  //                      hessian_pdot_p_FF    | hessian_pdot_pdot_FF    ]                   0  |  0 ]

  int i, j; // indexes used for number of agents (i is also used in other contexts)
  int k, l; // indexes used for the dimension of the space

  double SS[INTS];		  						// stores all sij computed
  double Del_p_[INTS][DIM]; 				// inter-agent displacement matrix
  double Dsig[INTS];    						// collection of sigma'
  double DDsig[INTS];    						// collection of sigma''

  // shape of the formation selected by the user
  //              		0      	  1				  2      	  3      	  4        	5
  /*double shape[] = {
		  								1.0,
		  								SQRT3,      1.0,
                      2.0, 	    SQRT3,  	  1.0,
		  								SQRT3, 	    2.0, 	    SQRT3, 	    1.0,
		  								1.0,      SQRT3,      2.0,      SQRT3,      1.0,
							 	 	 };*/

  /*double shape[] = {
		  								1.0,
		  								1.0,      1.0,
                      1.0, 	    1.0,  	  1.0,
		  								1.0, 	    1.0, 	    1.0, 	    1.0,
		  								1.0,      1.0,      1.0,      1.0,      1.0,
							 	 	 };*/

  // rectangular triangle in 2D
  /*double shape[] = {
                      1.0,
                      3.0/5.0,  4.0/5.0,
                   };*/

  // square in 3D
  /*double shape[] = {
                      1.0,
                      SQRT2,    1.0,
                      1.0,      SQRT2,    1.0,
                   };*/
  // pentagon in 2D
  //              		0      	  1				  2      	  3
  /*double shape[] = {
		  								1.0,
		  								DIAGP,    1.0,
                      DIAGP, 	  -1.0,  	  1.0,
		  								1.0, 	    DIAGP, 	  -1.0, 	  1.0,
							 	 	 };*/
	/*double shape[] = {
		  								1.0,
		  								DIAGP,    1.0,
                      DIAGP, 	  DIAGP,  	1.0,
		  								1.0, 	    DIAGP, 	  DIAGP, 	  1.0,
							 	 	 };*/
  // cube in 3D
  //                0       	  1				  2      	  3      	  4      	  5      	  6
  double shape[] = {
                    1.0,
                    -SQRT2,     1.0,
                    1.0,        SQRT2,    1.0,
                    1.0,        -SQRT2,   SQRT3,   -SQRT2,
                    SQRT2,      1.0,      SQRT2,   -SQRT3,    1.0,
                    -SQRT3,     -SQRT2,    1.0,      SQRT2,   SQRT2,    1.0,
                    SQRT2,      SQRT3,    -SQRT2,   1.0,      1.0,      -SQRT2,    1.0,
                 };
  /*double shape[] = {
                    1.0,
                    SQRT2,      1.0,
                    1.0,        SQRT2,    1.0,
                    1.0,        SQRT2,    SQRT3,    SQRT2,
                    SQRT2,      1.0,      SQRT2,    SQRT3,    1.0,
                    SQRT3,      SQRT2,    1.0,      SQRT2,    SQRT2,    1.0,
                    SQRT2,      SQRT3,    SQRT2,    1.0,      1.0,      SQRT2,    1.0,
                 };*/


	// computing FF(p) and inter-agent displacements
	double FF = 0.0;    			// formation function FF(sigma(||pi-pj||^2))
  if(do_l || do_a || do_Q || do_Qsafe){

    int ij = 0;							// index for interactions between different agents
    double del_p_ijk;  		  // projection on the k axis of the inter-agent displacement ij
    double sij;         		// squared inter-agent distance ij

    for(i = 1; i < nAg; i++){
      for(j = 0; j < i; j++){
      	// Just in case we would weight inter-agent distances
		    if(shape[ij] > 0.0){
			    sij = 0.0;
			    for(k = 0; k < DIM; k++){
			      del_p_ijk = x[i*DIM+k]-x[j*DIM+k];
			      Del_p_[ij][k] = del_p_ijk;
			      sij += del_p_ijk*del_p_ijk;
			    }
			    SS[ij] = sij;
		    	//    sigma(double s,       double d,   double *sig_s,   double *sig_s_s)
		    	FF += sigma(     sij,   shape[ij]*dd,       &Dsig[ij],        &DDsig[ij]);
				}
				else{
					Dsig[ij] = 0.0;
					DDsig[ij] = 0.0;
				}
				//FF += sigma(     sij,    dd,       &Dsig[ij],        &DDsig[ij]);
		    ij++;
      }
    }

  }

  // computing AA(dp) and inter-agent displacements
	double AA = 0.0;    			// formation function AA(||dpi-dpj||^2)
	double AlAl[INTS];		  				  // stores all al_ij computed
  double Del_dp_[INTS][DIM]; 				// inter-agent nonalignment matrix
  if(do_l || do_a || do_Q || do_Qsafe){

    int ij = 0;							// index for interactions between different agents
    double del_dp_ijk;  		  // projection on the k axis of the inter-agent displacement ij
    double al_ij;         		// squared inter-agent distance ij

    for(i = 1; i < nAg; i++){
      for(j = 0; j < i; j++){
      	// Just in case we would weight inter-agent velocities
		    if(shape[ij] > 0.0){
			    al_ij = 0.0;
			    for(k = 0; k < DIM; k++){
			      del_dp_ijk = x[NI+i*DIM+k]-x[NI+j*DIM+k];
			      Del_dp_[ij][k] = del_dp_ijk;
			      al_ij += del_dp_ijk*del_dp_ijk;
			    }
			    AlAl[ij] = al_ij;
		    	AA += 1.0*al_ij;   // Aij=1.0 is not yet programmed... it represents a weighting matrix for each link in the graph
				}
				else{
					AlAl[ij] = 0.0;
				}
				//AA += al_ij;
		    ij++;
      }
    }

  }


  // computing the gradient of FF(p) wrt the positions
  double nabla_p_FF[nAg][DIM];  // gradient of FF
  if(do_a){

  	int ij = 0;						      // index for interactions between different agents
  	double nabla_value_ijk;		  // value assigned to the gradient at each iteration

		// zero out nabla_p_FF
  	for(i = 0; i < nAg; i++){
  	  for(k = 0; k < DIM; k++){
	      nabla_p_FF[i][k] = 0.0;
	  	}
		}

  	// gradient is computed
  	for(i = 1; i < nAg; i++){
  		for(j = 0; j < i; j++){
  			if(shape[ij] > 0.0){
  				for(k = 0; k < DIM; k++){
  					nabla_value_ijk = Dsig[ij]*Del_p_[ij][k];
  					nabla_p_FF[i][k] += nabla_value_ijk;
  					nabla_p_FF[j][k] -= nabla_value_ijk;
					}
				}
				ij++;
			}
		}

  }


  // computing the gradient of AA(dp) wrt the positions
  double nabla_dp_AA[nAg][DIM];  // gradient of AA
  if(do_a){

  	int ij = 0;						      // index for interactions between different agents
  	double nabla_value_ijk;		  // value assigned to the gradient at each iteration

		// zero out nabla_dp_AA
  	for(i = 0; i < nAg; i++){
  	  for(k = 0; k < DIM; k++){
	      nabla_dp_AA[i][k] = 0.0;
	  	}
		}

  	// gradient is computed
  	for(i = 1; i < nAg; i++){
  		for(j = 0; j < i; j++){
  			if(shape[ij] > 0.0){
  				for(k = 0; k < DIM; k++){
  					nabla_value_ijk = Del_dp_[ij][k];
  					nabla_dp_AA[i][k] += nabla_value_ijk;
  					nabla_dp_AA[j][k] -= nabla_value_ijk;
					}
				}
				ij++;
			}
		}

  }


  // computing the Hessian of FF wrt the positions
  double H_p_p_FF[INTS][(DIM*(DIM+1))/2]; 	  // Hessian of FF
  if(do_Q || do_Qsafe){

  	int ij = 0;  // index for interactions between different agents
  	int kl;      // index for dimension of the space (used in matrix blocks)

  	// only off-diagonal blocks are computed
  	for(i = 1; i < nAg; i++){
  		for(j = 0; j < i; j++){
  			if(shape[ij] > 0.0){
  				kl = 0;
					for(k = 0; k < DIM; k++){
						for(l = 0; l <= k; l++){
							H_p_p_FF[ij][kl]  = 2.0*DDsig[ij]*(Del_p_[ij][k])*(-Del_p_[ij][l]);
							// prevent the matrix H be non-positive (semi)definite
	  					if(!do_Qsafe || Dsig[ij] > 0.0){
	  	  	  	  H_p_p_FF[ij][kl]  += Dsig[ij]*(k == l ? -1.0 : 0.0);
							}
							kl++;
						}
					}
				}
				ij++;
			}
		}

	}

	// computing the Hessian of AA wrt the velocities
  double H_dp_dp_AA[INTS][(DIM*(DIM+1))/2]; 	  // Hessian of AA
  if(do_Q || do_Qsafe){

  	int ij = 0;  // index for interactions between different agents
  	int kl;      // index for dimension of the space (used in matrix blocks)

  	// only off-diagonal blocks are computed
  	for(i = 1; i < nAg; i++){
  		for(j = 0; j < i; j++){
  			if(shape[ij] > 0.0){
  				kl = 0;
					for(k = 0; k < DIM; k++){
						for(l = 0; l <= k; l++){
							H_dp_dp_AA[ij][kl] = (k == l ? -1.0 : 0.0);
							kl++;
						}
					}
				}
				ij++;
			}
		}

	}


  // computing xB-xB_des
  if(do_l || do_a){
    for(k = 0; k < DIM; k++){
      dxB[k] = -xB_des[k];
      dxB[DIM+k] = -xB_des[DIM+k];
		}
		for(i = 0; i < nAg; i++){
			for(k = 0; k < DIM; k++){
				dxB[k] += x[i*DIM+k]/nAg;
				dxB[DIM+k] += x[(i+nAg)*DIM+k]/nAg;
			}
		}
  }



  if (do_a) {
    // zero out a
    for (i=0; i<NX; i++) {
      lxu_x_[i] = 0.0;
    }
  }

  if (do_b) {
    // zero out b
    for (i=0; i<NU; i++) {
      lxu_u_[i] = 0.0;
    }
  }

  if (do_Q || do_Qsafe) {
    // zero out Q = lxu_x_x
    for (i=0; i<NX*NX; i++) {
      lxu_x_x_[i] = 0.0;
    }
  }

  if (do_S) {
    // zero out S = lxu_x_u
    for (i=0; i<NX*NU; i++) {
      lxu_x_u_[i] = 0.0;
    }
  }

  if (do_R) {
    // zero out R = lxu_u_u
    for (i=0; i<NU*NU; i++) {
      lxu_u_u_[i] = 0.0;
    }
  }

  // l = nAg*|| xB - xB_des(t) ||^2_QB / 2 + || u ||^2_RD / 2 + FF(p) + AA(dp)
  if(do_l){
    l_fo = kF*FF/2.0 + kA*AA/2.0;
    l_st = 0.0;
    for(k = 0; k < DIM; k++) {
      l_st += nAg*q_p*dxB[k]*dxB[k]/2.0;
      l_st += nAg*q_v*dxB[DIM+k]*dxB[DIM+k]/2.0;
    }
    l_in = 0.0;
    for (i = 0; i < NU; i++){
      l_in += r_a*u[i]*u[i]/2.0;
    }
    lxu[0] = l_st+l_in+l_fo;
  }

  // a  = nAg*C'*QB*(xB-xB_des) + nabla_x_FF + nabla_x_AA
  if(do_a){
    for (i = 0; i < nAg; i++){
    	for(k = 0; k < DIM; k++){
    		lxu_x_[DIM*i+k] = q_p*dxB[k] + kF*nabla_p_FF[i][k];
    		lxu_x_[DIM*(nAg+i)+k] = q_v*dxB[DIM+k] + kA*nabla_dp_AA[i][k];
			}
    }
  }

  // b  = RD*u
  if(do_b){
    for(i = 0; i < nAg; i++){
      for(k = 0; k < DIM; k++){
        lxu_u_[i*DIM+k] = r_a*u[i*DIM+k];
      }
    }
  }

  // Q = nAg*C'*QB*C + H_x_x_FF + H_x_x_AA
  if (do_Q || do_Qsafe) {

    int ij = 0;          // index for interactions between different agents
  	int kl;              // index for dimension of the space (used in matrix blocks)
  	double Q_value_ijkl; // values of the Hessian to be assigned at each iteration (positions)
  	double dQ_value_ijkl; // values of the Hessian to be assigned at each iteration (velocities)

  	// nAg*C'*QB*C term (diagonal blocks)
    for(i = 0; i < nAg; i++){
      for(k = 0; k < DIM; k++){
        lxu_x_x(1+i*DIM+k,1+i*DIM+k) += q_p/nAg;
        lxu_x_x(1+(i+nAg)*DIM+k,1+(i+nAg)*DIM+k) += q_v/nAg;
      }
    }

    for(i = 0; i < nAg; i++){
      for(j = 0; j < i; j++){
        kl = 0;
        for(k = 0; k < DIM; k++){
          for(l = 0; l <= k; l++){

            // nAg*C'*QB*C term (off-diagonal blocks)
            if(k == l){
              lxu_x_x(1+i*DIM+k,1+j*DIM+l) += q_p/nAg;
              lxu_x_x(1+(i+nAg)*DIM+k,1+(j+nAg)*DIM+l) += q_v/nAg;
              // symmetry of the Hessian is implemented
              lxu_x_x(1+j*DIM+l,1+i*DIM+k) += q_p/nAg;
              lxu_x_x(1+(j+nAg)*DIM+l,1+(i+nAg)*DIM+k) += q_v/nAg;
            }

            // H_x_x_FF and H_x_x_AA terms
            if(shape[ij] > 0.0){
	            Q_value_ijkl = kF*H_p_p_FF[ij][kl];
							dQ_value_ijkl = kA*H_dp_dp_AA[ij][kl];
	            // 1a) value is added to block ij
	            lxu_x_x(1+i*DIM+k,1+j*DIM+l) += Q_value_ijkl;
	            lxu_x_x(1+nAg*DIM+i*DIM+k,1+nAg*DIM+j*DIM+l) += dQ_value_ijkl;
	            // 2a) value is added to block ji
	            lxu_x_x(1+j*DIM+l,1+i*DIM+k) += Q_value_ijkl;
	            lxu_x_x(1+nAg*DIM+j*DIM+l,1+nAg*DIM+i*DIM+k) += dQ_value_ijkl;
	            // 1b) value is subtracted from block ii
	            lxu_x_x(1+i*DIM+k,1+i*DIM+l) -= Q_value_ijkl;
	            lxu_x_x(1+nAg*DIM+i*DIM+k,1+nAg*DIM+i*DIM+l) -= dQ_value_ijkl;
	            // 2b) value is subtracted from block jj
	            lxu_x_x(1+j*DIM+k,1+j*DIM+l) -= Q_value_ijkl;
	            lxu_x_x(1+nAg*DIM+j*DIM+k,1+nAg*DIM+j*DIM+l) -= dQ_value_ijkl;
	            if(k != l){
	              // 3a) block symmetry is implemented in block ij
	              lxu_x_x(1+i*DIM+l,1+j*DIM+k) += Q_value_ijkl;
	              lxu_x_x(1+nAg*DIM+i*DIM+l,1+nAg*DIM+j*DIM+k) += dQ_value_ijkl;
	              // 4a) block symmetry is implemented in block ji
	              lxu_x_x(1+j*DIM+k,1+i*DIM+l) += Q_value_ijkl;
	              lxu_x_x(1+nAg*DIM+j*DIM+k,1+nAg*DIM+i*DIM+l) += dQ_value_ijkl;
	              // 3b) block symmetry is implemented in block ii
	              lxu_x_x(1+i*DIM+l,1+i*DIM+k) -= Q_value_ijkl;
	              lxu_x_x(1+nAg*DIM+i*DIM+l,1+nAg*DIM+i*DIM+k) -= dQ_value_ijkl;
	              // 4b) block symmetry is implemented in block jj
	              lxu_x_x(1+j*DIM+l,1+j*DIM+k) -= Q_value_ijkl;
	              lxu_x_x(1+nAg*DIM+j*DIM+l,1+nAg*DIM+j*DIM+k) -= dQ_value_ijkl;
	            }
						}
            kl++;
          }
        }
        ij++;
      }
    }

    /*printf("\n\nQ=[...\n");
    for(i = 1; i <= NS; i++){
      for(j = 1; j <= NS; j++){
        printf("%g ",lxu_x_x(i,j));
      }
      if(i == NS){
        printf("]");
      }
      printf(";\n");
    }*/

  }

  // S = 0
  if(do_S){
    // no lxu_x_u_ terms
  }

  // R = RD	(diagonal)
  if(do_R){
    for(i = 0; i < nAg; i++){
      for(k = 0; k < DIM; k++){
        lxu_u_u(1+i*DIM+k,1+i*DIM+k) = r_a;
      }
    }
  }

#undef lxu_x_x
#undef lxu_x_u
#undef lxu_u_u

#undef NX
#undef NU

}

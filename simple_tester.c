#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "sys_sizes.h"
#include "cost_params.h"
#include "sigma.c"

void main (){

  #define NX  (NS)
  #define NU  (NI)

  double lxu[1];
  double lxu_x_[NX], lxu_u_[NU];
  double lxu_x_x_[NX*NX], lxu_x_u_[NX*NU], lxu_u_u_[NU*NU];
	
  #define lxu_x_x(i,j)      lxu_x_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
  #define lxu_x_u(i,j)      lxu_x_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
  #define lxu_u_u(i,j)      lxu_u_u_[ ((i)-1) + ((j)-1)*(NU) ]  // columnwise
	
  #define dd     (1.0)   // desired distance to keep in formation
  #define dd2    (dd*dd)			        
  #define nAg2   (nAg*nAg)

  //double *x_des = wct;
  double x_des[] = {0.0, 0.0, 0.0,   0.0, 0.0, 0.0,   0.0, 0.0, 0.0,   1.0, 0.0, 0.0,   1.0, 0.0, 0.0,   1.0, 0.0, 0.0};
  double x[] =     {0.0, 0.0, 0.0,  -2.0, 1.0, 0.0,   0.0, -1.0, 0.0,  1.0, 0.0, 0.0,   1.0, 0.0, 0.0,   1.0, 0.0, 0.0};
  double u[] = {0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0 };
  double dx[NX], du[NU], dxB[2*DIM];
	
  bool do_l, do_a, do_b, do_Q, do_S, do_R, do_Qsafe;
  do_l = true;
  do_a = true;
  do_Q = true;
	
  double ll;
	
  // a  = (C'*QB*C+Qx)*(x-xd) + nabla_x_FF  | NS
  // b  = RD*u   						    | NI
  // Q = C'*QB*C+Qx + H_x_x_FF 		      | NSxNS
  // R = RD								  | NIxNI
  // S = 0 								  | NSxNI

  // computing FF(x) and its gradient and Hessian
  // hessian_x_x_FF = [hessian_p_p_FF       | hessian_p_pdot_FF          = [ hessian_p_p_FF  |  0
  //                   hessian_pdot_p_FF    | hessian_pdot_pdot_FF    ]                   0  |  0 ] 

  int i, j, k, l;
  
#define INTS ((nAg*(nAg-1))/2)			// number of interactions between different agents
#define S_INTS ((nAg*(nAg+1))/2)		// total number of possible interactions between agents
#define S_DIMS (DIM*(DIM+1)/2)			// total number of possible interactions between dimensions
  
  double SS[INTS];		  						// stores all sij computed
  double Del_p_[INTS][DIM]; 				// inter-agent displacement matrix
  double Dsig[INTS];    						// collection of sigma'
  double DDsig[INTS];    						// collection of sigma''

  double FF = 0.0;    			   			// formation function FF(sigma(||pi-pj||^2))
  double nabla_p_FF[nAg][DIM];     	// gradient of FF
  double H_p_p_FF[S_INTS][S_DIMS]; 	// Hessian of FF
  // bool do_Q_safe = false;
  // double H_p_p_FF_safe[S_INTS][S_DIMS];  // safe version for Hessian of FF



	// computing FF(x) and inter-agent displacements
  if (do_l || do_a || do_Q) {
	
    int ij = -1;						// index for interactions between different agents
    double del_p_ij_k;  		// projection on the k axis of the inter-agent displacement ij
    double sij;         		// squared inter-agent distance ij
    
    for(i = 1; i < nAg; i++){
      for(j = 0; j < i; j++){
	    ij++;
	    sij = 0.0;
	    for(k = 0; k < DIM; k++){
	      del_p_ij_k = x[i*DIM+k]-x[j*DIM+k];
	      Del_p_[ij][k] = del_p_ij_k;
	      sij += del_p_ij_k*del_p_ij_k; 
	    }
	    SS[ij] = sij;
			//    sigma(double s, double d, double *sig_s, double *sig_s_s)
	    FF += sigma(     sij,       dd,     &Dsig[ij],      &DDsig[ij]);
      }
    }
    
  }
  
    
  // computing the gradient of FF(x) wrt the positions
  if(do_a){
  	
  	int ij;								// index for interactions between different agents
  	double sign_Del_p_ij;	// sign of the ij-th displacement vector
  	
  	for(i = 0; i < nAg; i++){
  	  for(k = 0; k < DIM; k++){
	      nabla_p_FF[i][k] = 0.0;
	  	}
		}
  	
  	for(i = 0; i < nAg; i++){
  	  for(j = 0; j < nAg; j++){
		  	if( i != j){
		  		ij = (j < i ? ((i*(i-1))/2+j) : ((j*(j-1))/2+i) );
		  		sign_Del_p_ij = (j < i ? 1.0 : -1.0);
		  		for(k = 0; k < DIM; k++){
		  			nabla_p_FF[i][k] += Dsig[ij]*(2.0*sign_Del_p_ij*Del_p_[ij][k]);
		  		}
				}
			}
  	}
  	
  }
  
  
  // computing the Hessian of FF(x) wrt the positions
  if(do_Q){
  	
  	int ij;									// index for interactions between different agents
  	int kl_block;						// index for interactions between dimensions for Hessian blocks
  	int ij_block, ii_block; // off-diagonal and on-diagonal indexes for Hessian blocks
  	
  	for(i = 0; i < nAg; i++){
  		
  		// initialization of the current on-diagonal block wrt the i-th row
			ii_block = (i*(i+3))/2;
			kl_block = -1;
  		for(k = 0; k < DIM; k++){
		  	for(l = 0; l <= k; l++){
		  		H_p_p_FF[ii_block][++kl_block] = 0.0;
		  	}
		  }
  		
  	  for(j = 0; j < nAg; j++){
		  	if( i != j){
		  		
		  		ij = (j < i ? ((i*(i-1))/2+j) : ((j*(j-1))/2+i) );
		  		ij_block = (j < i ? ij+i : ij+j);
		  		kl_block = -1;
		  		
		  		// off-diagonal block computation
		  		for(k = 0; k < DIM; k++){
		  			for(l = 0; l <= k; l++){
		  				kl_block++;
		  				// block and matrix symmetries are imposed by construction
							if(j > i){
								H_p_p_FF[ij_block][kl_block] = DDsig[ij]*(2.0*Del_p_[ij][k])*(-2.0*Del_p_[ij][l]);
								// prevent the matrix be non-positive (semi)definite
		  					if(Dsig[ij] >= 0.0){
	    	  	  		H_p_p_FF[ij_block][kl_block] += Dsig[ij]*(k == l ? -2.0 : 0.0);
								}
							}
							// blocks on the diagonal are computed as minus the summation of the off-diagonal terms wrt its row					
		  				H_p_p_FF[ii_block][kl_block] -= H_p_p_FF[ij_block][kl_block];
						}
		  		}
		  		
				}
			}
  	}
  	
	}
    
    // I haven't thought yet about a nice theory to perform this task, I'm just using heuristics.
    if(do_a){
      double rescale = 100000.0;
  	  double min_dist = dd/rescale;
  	  double dt_ = 0.01;
  	  int ij;
  	  for(ij = 0; ij < ((nAg-1)*nAg)/2; ij++){
  	  	if(SS[ij] < min_dist){
	  	  for(i = 0; i < nAg; i++){
	  	    for(k = 0; k < DIM; k++){
	  	  	  nabla_p_FF[i][k] += min_dist/dt_*rand()/RAND_MAX;
		    }
		  }
		  break;
	    }
	  }
	  
	}
	
  
    
    for(i = 0; i < nAg; i++){
      for(j = 0; j < DIM; j++)
        printf("%f ", nabla_p_FF[i][j]);
      printf("\n");
	}
	printf("\n________________________________________________________________________\n");
	
	
	for (i=1; i<=(nAg*(nAg+1))/2; i++){
    for(j=1; j<=(DIM*(DIM+1))/2; j++){
    	printf("%f ", H_p_p_FF[i][j]);
    }
    printf("\n");
  }
  printf("\n________________________________________________________________________\n");
	
	int ih, jh, kh, lh;
    for (i=1; i<=NU; i++){
      for(j=1; j<=NU; j++){
	    	kh = (i-1)%DIM;
	    	lh = (j-1)%DIM;
	      ih = (i-1)/DIM;
	      jh = (j-1)/DIM;
	      if(i > j)
	        ih = (ih*(ih+1))/2;
	      else
	        jh = (jh*(jh+1))/2;
	      if(kh >= lh)
	      	kh = (kh*(kh+1))/2;
	      else
	      	lh = (lh*(lh+1))/2;
	      printf("%f ", H_p_p_FF[ih+jh][kh+lh]);
      }
      printf("\n");
	}
	printf("\n________________________________________________________________________\n");
	
	
	for(i = 0; i < nAg*(nAg-1)/2; i++){
	  printf("%f ", DDsig[i]);
	}
	printf("\n________________________________________________________________________\n");
    
    
    for(i = 0; i < nAg*(nAg-1)/2; i++){
      for(j = 0; j < DIM; j++)
	    printf("%f ", Del_p_[i][j]);
	  printf("\n");
	}
	printf("\n");
	printf("\n________________________________________________________________________\n");
	
	for(i = 0; i < nAg*(nAg-1)/2; i++){
	  printf("%f ", SS[i]);
	}
	printf("\n________________________________________________________________________\n");
	
  
    

  // make states and controls more accessible
  //   not needed here


 
    /*for(i = 0; i<NX; i++)
      dx[i] = x[i]-x_des[i];
    for(i = 0; i<2*DIM; i++)
      dxB[i] = 0.0;
    for(i=0; i<NU; i++){
      j = i%DIM;
      dxB[j] += dx[i]/nAg;
      dxB[DIM+j] += dx[NU+i]/nAg;
    }
  

 
    for (i=0; i<NX; i++) {
      lxu_x_[i] = 0.0;
    }
  

 
    // zero out b
    for (i=0; i<NU; i++) {
      lxu_u_[i] = 0.0;
    }
  


    // zero out Q = lxu_x_x
    for (i=0; i<NX*NX; i++) {
      lxu_x_x_[i] = 0.0;
    }
  

 
    // zero out S = lxu_x_u
    for (i=0; i<NX*NU; i++) {
      lxu_x_u_[i] = 0.0;
    }
  

 
    // zero out R = lxu_u_u
    for (i=0; i<NU*NU; i++) {
      lxu_u_u_[i] = 0.0;
    }
  

  // l = || x - x_des(t) ||^2_C'*QB*C / 2 + || u ||^2_R / 2 + FF(x)
  
    ll = FF;
    for(i=0; i<2*DIM; i++)
      ll += qB[i]*dxB[i]*dxB[i]/2.0;
    for (i=0; i<NU; i++)
      ll += RD[i]*u[i]*u[i]/2.0;
    lxu[0] = ll;
  
    printf("%f ", lxu[0]);
    printf("\n___________________________________________________________________________\n");

  // a  = C'*QB*C*(x-x_des) + nabla_x_FF  
 
    for (i=0; i<NX; i++){
      j = i%DIM;
      if(i < NU)
	    lxu_x_[i] = qB[j]/nAg2*dx[i]+nabla_p_FF[i/DIM][j];
      else
	    lxu_x_[i] = qB[DIM+j]/nAg2*dx[NU+i];
	  printf("%f ", lxu_x_[i]);
    }	
    
    printf("\n___________________________________________________________________________\n");

  // b  = RD*u   						      
  
    for (i=0; i<NU; i++)
      lxu_u_[i] = RD[i]*u[i];
  

  // Q = C'*QB*C + H_x_x_FF  
 
    for (i=1; i<=NX; i++){
      for(j=1; j<=NX; j++){
	    kh = (i-1)%DIM;
	    lh = (j-1)%DIM;
	    if(i <= NU && j <= NU){
	      ih = (i-1)/DIM;
	      jh = (j-1)/DIM;
	      if(i > j)
	        ih = ih*(ih+1)/2;
	      else
	        jh = jh*(jh+1)/2;
	      lxu_x_x(i,j) = H_p_p_FF[ih+jh][kh*(kh+1)/2+lh];
	      lxu_x_x(i,j) += (kh == lh ? qB[kh]/nAg2 : 0.0);
	    }
	    else if (i > NU && j > NU && kh == lh) {
	      lxu_x_x(i,j) = qB[kh+DIM]/nAg2;
	    }
	    printf("%f ", lxu_x_x(i,j));
      }
      printf("\n");
    }

  // S = 0 
  
    // no lxu_x_u_ terms 


  // R = RD	(diagonal)

    for (i=1; i<=NU; i++)
      lxu_u_u(i,i) = RD[i-1];*/
  

}

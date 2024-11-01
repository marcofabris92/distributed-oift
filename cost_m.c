/*
 *
 * cost_m.c
 *
 *   .MEX function to evaluate
 *
 *      system dynamics and derivatives
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt);
 *
 * to evaluate, we call
 *
  cost(x,u,wlt, ders, lxu,    lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
                ders: lxu(1), a(2),  b(4),   Q(8),    S(16),   R(32)
 *
 *  ders is determined from nargout
 *
 * JH, Boulder, apr 12
 * checked oct 13
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#include <stdio.h>
#include <math.h>
// #include <memory.h>
#include "mex.h"

#include "sys_sizes.h"
#include "cost.c"

char sys_name[] = "cost";

#define  NX    (NS)
#define  NU    (NI)

/* Input Arguments */

#define	XX        prhs[0]
#define	UU        prhs[1]
#define WLT       prhs[2]
#define do_Qsafe  prhs[3]

/* Output Arguments */

#define	LXU        plhs[0]
#define LXU_X      plhs[1]
#define LXU_U      plhs[2]
#define LXU_X_X    plhs[3]
#define LXU_X_U    plhs[4]
#define LXU_U_U    plhs[5]

#if !defined(max)
#define	max(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(min)
#define	min(A, B)	((A) < (B) ? (A) : (B))
#endif

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double *x, *u, *wlt, *lxu;
  double
    *lxu_x_ = NULL, *lxu_u_ = NULL,
    *lxu_x_x_ = NULL, *lxu_x_u_ = NULL, *lxu_u_u_ = NULL;
  int ders, nx, nu, nwlt;
  unsigned int	m, n;

  char errMsg[256];

  /*
   * things to check:
   *
   *  nrhs > 1 always        (need x and u)
   *  nrhs > 2 if NWL > 0    (need wlt)
   *  sizes of x,u,wlt are correct
   *
   */

  // printf("nrhs: %d, nlhs: %d\n",nrhs,nlhs);

  // nrhs > 1 always        (need x and u)
  if (nrhs < 2) {
    sprintf(errMsg, "%s_m requires x and u:\n  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = %s_m(x, u, wlt);", sys_name, sys_name);
    mexErrMsgTxt(errMsg);
  }

  // nrhs > 2 if NWL > 0    (need wlt)
  if ( NWL > 0 && nrhs < 3 ) {
    sprintf(errMsg, "%s_m requires a wlt input of size %d:\n  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = %s_m(x, u, wlt);", sys_name, NWL, sys_name);
    mexErrMsgTxt(errMsg);
  }

  // get input matrices/vectors

  x   = mxGetPr(XX);
  u   = mxGetPr(UU);

  // *** ADD *** size check on x and u

  if (nrhs > 2) {
    nwlt = mxGetNumberOfElements(WLT);
    if (nwlt != NWL) {
      sprintf(errMsg, "%s_m: wlt has %d elements but should have %d:\n  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = %s_m(x, u, wlt);", sys_name, nwlt, NWL, sys_name);
      mexErrMsgTxt(errMsg);
    }
    wlt = mxGetPr(WLT);
  }

  // determine number of derivs wanted
  ders = 1;                           // lxu
  if (nlhs > 1) ders += 2 + 4;        // A & B
  if(nlhs > 3 && !do_Qsafe){
    ders += 8 + 16 + 32;  // Q & S & R
  }
  else{
    ders += 16 + 32 + 64;  // Qsafe & S & R
  }

  // create lxu output
  LXU = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  lxu = mxGetPr(LXU);

  if (nlhs > 1) {
    // create first derivative matrices for output
    LXU_X = mxCreateDoubleMatrix((mwSize)1, (mwSize)NX, mxREAL);
    LXU_U = mxCreateDoubleMatrix((mwSize)1, (mwSize)NU, mxREAL);

    lxu_x_ = mxGetPr(LXU_X);
    lxu_u_ = mxGetPr(LXU_U);
  }

  if (nlhs > 3) {
    // create second derivative matrices for output
    LXU_X_X = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NX, mxREAL);
    LXU_X_U = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NU, mxREAL);
    LXU_U_U = mxCreateDoubleMatrix((mwSize)NU, (mwSize)NU, mxREAL);

    lxu_x_x_ = mxGetPr(LXU_X_X);
    lxu_x_u_ = mxGetPr(LXU_X_U);
    lxu_u_u_ = mxGetPr(LXU_U_U);
  }

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  cost(x,u,wlt,ders,lxu, lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);

#undef NX
#undef NU
}

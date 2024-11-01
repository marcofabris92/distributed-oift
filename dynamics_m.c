/*
 *
 * dynamics_m.c
 *
 *   .MEX function to evaluate
 *
 *      system dynamics and derivatives
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
 [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = dynamics_m(x, u, wt, q);
 *
 * to evaluate, we call
 *
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
         ders:     dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
 *
 *  ders is determined from nargout
 *
 *  x and u are always required
 *  wt is required if NW > 0
 *  q  is required if nlhs > 4
 *
 *  sizes will be checked
 *
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#include <stdio.h>
#include <math.h>
// #include <memory.h>
#include "mex.h"

#include "sys_sizes.h"
#include "dynamics.c"

// "sys_sizes.h" defines
// NS   - number of system states
// NI   - number of system inputs
// NW   - number of exogenous inputs for dynamics
// NWL  - number of exogenous inputs for cost
// NWC  - number of exogenous inputs for constraints

char sys_name[] = "dynamics";

#define  NX    (NS)
#define  NU    (NI)

/* Input Arguments */

#define	XX     prhs[0]
#define	UU     prhs[1]
#define WT     prhs[2]
#define	QQ     prhs[3]

/* Output Arguments */

#define	DX         plhs[0]
#define YY         plhs[1]
#define FXU_X      plhs[2]
#define FXU_U      plhs[3]
#define Q_FXU_X_X  plhs[4]
#define Q_FXU_X_U  plhs[5]
#define Q_FXU_U_U  plhs[6]

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
  double *x, *u, *wt = NULL, *dx, *y = NULL;
  double
    *fxu_x_ = NULL, *fxu_u_ = NULL, 
    *q = NULL, *q_fxu_x_x_ = NULL, *q_fxu_x_u_ = NULL, *q_fxu_u_u_ = NULL;
  int ders, der, nx, nu, nwt, nq;
  unsigned int	m, n;
  int i;

  char errMsg[256];

  /*
   * things to check:
   *
   *  nrhs > 1 always        (need x and u)
   *  nrhs > 2 if NW > 0     (need wt)
   *  nrhs > 3 if nlhs > 4   (need q for 2nd derivs)
   *  sizes of x,u,wt,q are correct
   *
   */

  // printf("nrhs: %d, nlhs: %d\n",nrhs,nlhs);

  // nrhs > 1 always        (need x and u)
  if (nrhs < 2) {
    sprintf(errMsg, "%s_m requires x and u:\n  [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = %s_m(x, u, wt, q);", sys_name, sys_name);
    mexErrMsgTxt(errMsg);
  }

  // nrhs > 2 if NW > 0    (need wt)
  if ( NW > 0 && nrhs < 3 ) {
    sprintf(errMsg, "%s_m requires a wt input of size %d:\n  [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = %s_m(x, u, wt, q);", sys_name, NW, sys_name);
    mexErrMsgTxt(errMsg);
  }

  // nrhs > 3 if nlhs > 4   (need q for 2nd derivs)
  if ( nlhs > 4 && nrhs < 4 ) {
    sprintf(errMsg, "%s_m: 2nd derivs require q:\n  [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = %s_m(x, u, wt, q);", sys_name, sys_name);
    mexErrMsgTxt(errMsg);
  }

  // get input matrices/vectors
  
  x  = mxGetPr(XX);
  u  = mxGetPr(UU);

  // *** ADD *** size check on x and u

  if (nrhs > 2) {
    nwt = mxGetNumberOfElements(WT);
    if (nwt != NW) {
      sprintf(errMsg, "%s_m: wt has %d elements but should have %d:\n  [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = %s_m(x, u, wt, q);", sys_name, nwt, NW, sys_name);
      mexErrMsgTxt(errMsg);
    }
    wt = mxGetPr(WT);
  }

  if (nrhs > 3) {
    // *** ADD *** size check on q
    q = mxGetPr(QQ);
  }

  // determine number of derivs wanted
  ders = 1;
  der  = 1;
  for (i=1; i<nlhs; i++) {
    der  *= 2;
    ders += der;
  }

  // create dx vector output
  DX = mxCreateDoubleMatrix((mwSize)NX, (mwSize)1, mxREAL);
  dx = mxGetPr(DX);

  if (nlhs > 1) {
    YY = mxCreateDoubleMatrix((mwSize)NO, (mwSize)1, mxREAL);
    y = mxGetPr(YY);
  }

  if (nlhs > 2) {
    FXU_X = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NX, mxREAL);
    fxu_x_ = mxGetPr(FXU_X);
  }

  if (nlhs > 3) {
    FXU_U = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NU, mxREAL);
    fxu_u_ = mxGetPr(FXU_U);
  }

  if (nlhs > 4) {
    Q_FXU_X_X = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NX, mxREAL);
    q_fxu_x_x_ = mxGetPr(Q_FXU_X_X);
  }

  if (nlhs > 5) {
    Q_FXU_X_U = mxCreateDoubleMatrix((mwSize)NX, (mwSize)NU, mxREAL);
    q_fxu_x_u_ = mxGetPr(Q_FXU_X_U);
  }

  if (nlhs > 6) {
    Q_FXU_U_U = mxCreateDoubleMatrix((mwSize)NU, (mwSize)NU, mxREAL);
    q_fxu_u_u_ = mxGetPr(Q_FXU_U_U);
  }
  

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  dynamics(x,u,wt,ders,dx,y, fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);

#undef NX
#undef NU
}

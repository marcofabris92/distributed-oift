/*
 *
 * sigma_m.c
 *
 *   .MEX function to evaluate
 *
 *      sigma(s) function and derivatives
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
  [sig, sig_s, sig_s_s] = sigma_m(s,d);
 *
 * to evaluate, we call
 *
 *  double sigma(double s, double d, double *sig_s, double *sig_s_s)
 *
 *
 * JH jul18 boulder
 */

#include <stdio.h>
#include <math.h>
// #include <memory.h>
#include "mex.h"

#include "sigma.c"

char sys_name[] = "sigma";

/* Input Arguments */

#define	SS     prhs[0]
#define	DD     prhs[1]

/* Output Arguments */

#define	SIG        plhs[0]
#define SIG_S      plhs[1]
#define SIG_S_S    plhs[2]

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
  double *s_, *d_, *sig_, *sig_s_, *sig_s_s_;

  // char errMsg[256];

  // get input matrices/vectors
  s_ = mxGetPr(SS);
  d_ = mxGetPr(DD);

  // create outputs
  SIG = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  sig_ = mxGetPr(SIG);

  SIG_S = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  sig_s_ = mxGetPr(SIG_S);

  SIG_S_S = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  sig_s_s_ = mxGetPr(SIG_S_S);

  *sig_ = sigma(*s_,*d_,sig_s_,sig_s_s_);
}

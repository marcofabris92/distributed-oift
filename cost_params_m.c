/*
 *
 * cost_params_m.c
 *
 *   .MEX function to obtain the
 *
 *     cost paramaters  
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
 *   [Qreg, Rreg, Q, R] = cost_params_m();
 *
 *
 */

#include "mex.h"

#include "sys_sizes.h"
#include "QR_params.h"
#include "cost_params.h"

/* "QR_params.h" defines */
/* Qr - regulator state cost */
/* Rr - regulator input cost */

/* "cost_params.h" defines */
/* QD - State Cost */
/* RD - Input Cost */

/* no input arguments */

/* Output Arguments */

#define	_Qr   plhs[0]
#define _Rr   plhs[1]
#define _QD   plhs[2]
#define _RD   plhs[3]

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double *qr, *rr, *qd, *rd;
  int i;


  // create and assign output variables
  _Qr = mxCreateDoubleMatrix((mwSize)NS, (mwSize)NS, mxREAL);
  qr = mxGetPr(_Qr);
  for(i=0;i<NS;i++)
  {
    *qr = Qr[i];
    qr += NS+1;
  }

  if (nlhs > 1) {
    _Rr = mxCreateDoubleMatrix((mwSize)NI, (mwSize)NI, mxREAL);
    rr = mxGetPr(_Rr);
    for(i=0;i<NI;i++)
    {
      *rr = Rr[i];
      rr += NI+1;
    }
  }

  if (nlhs > 2) {
    _QD = mxCreateDoubleMatrix((mwSize)NS, (mwSize)NS, mxREAL);
    qd = mxGetPr(_QD);
    for(i=0;i<NS;i++)
    {
      *qd = QD[i];
      qd += NS+1;
    }
  }

  if (nlhs > 3) {
    _RD = mxCreateDoubleMatrix((mwSize)NI, (mwSize)NI, mxREAL);
    rd = mxGetPr(_RD);
    for(i=0;i<NI;i++)
    {
      *rd = RD[i];
      rd += NI+1;
    }
  }
}


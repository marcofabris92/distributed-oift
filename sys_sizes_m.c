/*
 *
 * sys_sizes_m.c
 *
 *   .MEX function to obtain the
 *
 *      system sizes
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
 *   [NS, NI, NO, NW, NWL, NWC] = sys_sizes_m();
 *
 *
 */

#include "mex.h"

#include "sys_sizes.h"

// "sys_sizes.h" defines
// NS   - number of system states
// NI   - number of system inputs
// NW   - number of exogenous inputs for dynamics
// NWL  - number of exogenous inputs for cost
// NWC  - number of exogenous inputs for constraints

/* no input arguments */

/* Output Arguments */

#define	_NS   plhs[0]
#define _NI   plhs[1]
#define _NO   plhs[2]
#define _NW   plhs[3]
#define _NWL  plhs[4]
#define _NWC  plhs[5]


void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double *ns, *ni, *no, *nw, *nwl, *nwc;


  // create and assign output variables
  _NS = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  ns = mxGetPr(_NS);
  *ns = NS;

  if (nlhs > 1) {
    _NI = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    ni = mxGetPr(_NI);
    *ni = NI;
  }

  if (nlhs > 2) {
    _NO = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    no = mxGetPr(_NO);
    *no = NO;
  }

  if (nlhs > 3) {
    _NW = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    nw = mxGetPr(_NW);
    *nw = NW;
  }

  if (nlhs > 4) {
    _NWL = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    nwl = mxGetPr(_NWL);
    *nwl = NWL;
  }

  if (nlhs > 5) {
    _NWC = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    nwc = mxGetPr(_NWC);
    *nwc = NWC;
  }

}

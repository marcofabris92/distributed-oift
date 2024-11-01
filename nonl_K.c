#define S_FUNCTION_NAME  nonl_K

/*
 * nonl_K.c
 *
 * projection operator + cost evaluation for
 *
 *   dynamics.c  system
 *
 * with
 *
 *   cost.c incremental cost
 *
 *  states
 *
 *    x
 *    cost
 *
 *  inputs
 *
 *    xi = [alf mu]     traj to track
 *    Kr   feedback gain for projection operator
 *    wt   exogenous input (if any) for dynamics
 *    wlt  exogenous input (if any) for cost
 *
 *  outputs
 *
 *    ux
 *
 *  parameters
 *
 *    (none)     [[ e.g., Rhg        hand of god weight ]]
 *
 * John Hauser
 * Nov/Dec 06
 * ...
 * JH, Boulder, Dec 11, Jan 12
 *
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#include <math.h>

// system description: dynamics with derivatives
#include "sys_sizes.h"
#include "dynamics.c"

// cost description, including derivatives
#include "cost.c"

// "sys_sizes.h" defines
// NS   - number of system states
// NI   - number of system inputs
// NO   - number of system outputs
// NW   - number of exogenous inputs for dynamics
// NWL  - number of exogenous inputs for cost

// number of S-function states, inputs, ...
#define N_STATES	(NS+1)
#define N_INPUTS	(NS+NI+NI*NS+NW+NWL)
#define N_OUTPUTS	(NI+NO)
#define N_PARAMS	(0)

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    /* See sfuntmpl.doc for more details on the macros below */

    ssSetNumSFcnParams(S, N_PARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, N_STATES);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, N_INPUTS);
    // ssSetInputPortDirectFeedThrough(S, 0, 0);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, N_OUTPUTS);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* set initial conditions with  simset  */
#if defined(MDL_INITIALIZE_CONDITIONS)
#undef MDL_INITIALIZE_CONDITIONS
#endif


// compute projection operator input
//   ux = mu(t) + K(t) [ alf(t) - x ]
static void fdbk(double *x, double *u, double *ux)
{

#define K(i,j)	K_[ ((i)-1)*NS + ((j)-1) ]	// rowwise

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS, *K_ = mu+NI;  // *wt = K_+NI*NS, *wlt = wt+NW;
  int i,j;

  /* compute the projection operator feedback */

  for (i=1; i<=NI; i++) {
    ux[i-1] = mu[i-1];
    for (j=1; j<=NS; j++) {
      ux[i-1] += K(i,j)*(alf[j-1] - x[j-1]);
    }
  }

#undef K
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
  real_T	*y	    = ssGetOutputPortRealSignal(S,0);
  real_T	*x	    = ssGetContStates(S);
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double *param0 = mxGetPr(ssGetSFcnParam(S,0));

  double u[N_INPUTS];
  int i;

  /* outputs include feedback control and system output(s) */
  double *ux = y, *yy = ux+NI;

  // pointers to 'input structure'
  //   we may need exogenous input w(t) to evaluate system output
  double *alf = u, *mu = alf+NS, *K_ = mu+NI, *wt = K_+NI*NS, *wlt = wt+NW;

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

  // ux = mu(t) + K(t) [ alf(t) - x ]    (( ux is the system inpt u ))
  fdbk(x, u, ux);

  // get system output if defined
  if (NO > 0) {
    dynamics(x,ux,wt, 2, NULL, yy,  NULL,NULL,  NULL, NULL,NULL,NULL);
  }

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

}

/* no model update */
#if defined(MDL_UPDATE)
#undef MDL_UPDATE
#endif

#define MDL_DERIVATIVES
static void mdlDerivatives(SimStruct *S)
{
  real_T            *dx     = ssGetdX(S);
  real_T            *x      = ssGetContStates(S);
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double *param0 = mxGetPr(ssGetSFcnParam(S,0));

  double u[N_INPUTS];
  int i, j, k;

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS, *K_ = mu+NI, *wt = K_+NI*NS, *wlt = wt+NW;

  double ux[NI];

  double lxu[0];  // incremental cost

  /* copy inputs to uu[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

  // ux = mu(t) + K(t) [ alf(t) - x ]    (( ux is the system inpt ))
  fdbk(x, u, ux);

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  dynamics(x,ux,wt, 1, dx, NULL,  NULL,NULL,  NULL, NULL,NULL,NULL);

  // short circuit for unbounded system
  dynamics_unbdd(x,ux,wt, dx);

  cost(x,ux,wlt, 1, lxu,  NULL,NULL,  NULL,NULL,NULL);

  dx[NS] = lxu[0];
}

static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

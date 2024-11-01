#define S_FUNCTION_NAME linK_zv

/*
 * linK_zv.c
 *
 * run linearized system
 *   for descent direction (or neighboring extremal)
 *   with 'optimal' feedback Ko 
 *   (determined by P, etc.)
 * for
 *
 *   dynamics.c  system
 *
 * with
 *
 *   cost.c incremental cost
 *
 *
 * John Hauser
 * Mar 04
 * ...
 * Dec 11
 *
 *  states
 *	z(NS)	linearized state 
 *	Dg(xi) . zeta                   lin approx to cost
 *	D^2g(xi) . (zeta, zeta)	        quad approx to cost
 *
 *  inputs
 *	xi(NS+NI)       alf(NS) mu(NI)  - trajectory for linearization
 *	Ko(NS*NI)       optimal feedback
 *	vo(NI)          optimal feedforward input
 *	q(NS)           costate traj (for 2nd order approx)
 *	wt(NW)         exogenous system input w(t)
 *	wlt(NWL)        exogenous cost input w_c(t) (e.g., x_des, u_des)
 *
 *  outputs
 *	v(NI)	linearized control
 *
 *  parameters
 *	none
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

#define NS2		(NS*NS)
#define NN12		((NS*(NS+1))/2)

#define N_STATES	(NS+1+1)
#define N_INPUTS	(NS+NI + NI*NS + NI + NS + NW + NWL)
#define N_OUTPUTS	(NI)
#define N_PARAMS	(0)


static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, N_PARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, N_STATES);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, N_INPUTS);
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

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* set initial conditions with  set options */
#if defined(MDL_INITIALIZE_CONDITIONS)
#undef MDL_INITIALIZE_CONDITIONS
#endif


/*
 * most matrices are stored COLUMNwise
 *  THE exception is the feedback gain matrix K
 *  which is stored ROWwise
 *
 * in each case, we set up macros to allow matlab-like indexing
 *   [ (1..n) vs (0..n-1) ]
 */

static void opt_v(double *v, double *z, double *u)
{
#define Ko(i,j)	Ko_[ ((i)-1)*NS + ((j)-1) ]	// rowwise

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS,
    *Ko_ = mu+NI, *vo = Ko_+(NI*NS), *q = vo+NI,
    *wt = q+NS, *wlt = wt+NW;

  int i, j, k;

  /* v = -Ko * z + vo  */
  for (i=1; i<=NI; i++) {
    v[i-1] = vo[i-1];
    for (j=1; j<=NS; j++) {
      v[i-1] -= Ko(i,j)*z[j-1];
    }
  }

#undef Ko
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
  real_T	*y	= ssGetOutputPortRealSignal(S,0);
  real_T	*x	= ssGetContStates(S);	
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double     *param0 = mxGetPr(ssGetSFcnParam(S,0));

  double u[N_INPUTS];
  int i;

  double *v=y, *z=x;

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

  /* output the control u and feedback Kux */
  opt_v(v, z, u);
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
  // double         *params = mxGetPr(ssGetSFcnParam(S,0));

  double u[N_INPUTS];
  int i, j, k;

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS,
    *Ko_ = mu+NI, *vo = Ko_+(NI*NS), *q = vo+NI,
    *wt = q+NS, *wlt = wt+NW;

  double *z = x;

  double A_[NS2], B_[NS*NI], v[NI];
  double a[NS], b[NI];
  double Dl_zv, QRzv2;

  double   lxu_x_x_[NS*NS],   lxu_x_u_[NS*NI],   lxu_u_u_[NI*NI];
  double q_fxu_x_x_[NS*NS], q_fxu_x_u_[NS*NI], q_fxu_u_u_[NI*NI];

  /* index A & P with 1:NS  (NOT 0:NS-1!) */
#define A(i,j)	A_[((i)-1) + NS*((j)-1)]
#define B(i,j)	B_[ ((i)-1) + ((j)-1)*(NS) ]	// columnwise

#define   lxu_x_x(i,j)      lxu_x_x_[((i)-1) + NS*((j)-1)]
#define   lxu_x_u(i,j)      lxu_x_u_[((i)-1) + NS*((j)-1)]
#define   lxu_u_u(i,j)      lxu_u_u_[((i)-1) + NI*((j)-1)]

#define q_fxu_x_x(i,j)    q_fxu_x_x_[((i)-1) + NS*((j)-1)]
#define q_fxu_x_u(i,j)    q_fxu_x_u_[((i)-1) + NS*((j)-1)]
#define q_fxu_u_u(i,j)    q_fxu_u_u_[((i)-1) + NI*((j)-1)]

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  // get derivs of vector field and cost about  (alf(.),mu(.))
  // need A B a b Q S R

  // a(2), b(4), Q(8), S(16), R(32)
  cost(    alf,mu,wlt,
           2+4+8+16+32,NULL, a, b,       lxu_x_x_,  lxu_x_u_,  lxu_u_u_);
  // A(4), B(8), Q(16), S(32), R(64)
  dynamics(alf,mu,wt,
           4+8+16+32+64,NULL,NULL, A_,B_, q, q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);

  /* zdot = A z + B v */

  opt_v(v, z, u);

  for (i=1; i<= NS; i++) {
    dx[i-1] = 0.0;
    for (j=1; j<=NS; j++) {
      dx[i-1] += A(i,j)*z[j-1];
    }
  }

  for (i=1; i<= NS; i++) {
    for (j=1; j<=NI; j++) {
      dx[i-1] += B(i,j)*v[j-1];
    }
  }

  /* Dg(xi) . zeta */
  Dl_zv = 0.0;
  for (i=1; i<=NS; i++) {
    Dl_zv += a[i-1]*z[i-1];
  }
  for (i=1; i<=NI; i++) {
    Dl_zv += b[i-1]*v[i-1];
  }
  dx[NS] = Dl_zv;


  /* D^2g(xi) . (zeta, zeta) */

  QRzv2 = 0.0;
  for (i=1; i<= NS; i++){
    for (j=1; j<= NS; j++){
      QRzv2 += z[i-1]*(lxu_x_x(i,j)+q_fxu_x_x(i,j))*z[j-1];
      }
  }

  for (i=1;i<=NS;i++) {
    for (j=1;j<=NI;j++) {
      QRzv2 += 2.0*z[i-1]*(lxu_x_u(i,j)+q_fxu_x_u(i,j))*v[j-1];
    }
  }

  for (i=1; i<=NI; i++) {
    for (j=1; j<=NI; j++){
      QRzv2 += v[i-1]*(lxu_u_u(i,j)+q_fxu_u_u(i,j))*v[j-1];
    }
  }

  dx[NS+1] = QRzv2;

#undef A
#undef B

#undef q_fxu_x_x
#undef q_fxu_x_u
#undef q_fxu_u_u

#undef lxu_x_x
#undef lxu_x_u
#undef lxu_u_u
}

/*
 * mdlTerminate - called when the simulation is terminated.
 *
 * In this function, you should perform any actions that are necessary
 * at the termination of a simulation.  For example, if memory was allocated
 * in mdlInitializeConditions, this is the place to free it.
 */

static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

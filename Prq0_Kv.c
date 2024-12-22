// #define SECOND_ORDER_DESCENT
// #define S_FUNCTION_NAME Prq2_Kv

/*
 * Prq0_Kv.c
 *
 * MASTER S-function to compute
 *   P, r, and q in the calculation of
 * BOTH
 *   FIRST   (without D^2 P terms)
 * and
 *   SECOND  (includes D^2 P terms)
 * ORDER descent directions
 *
 * for
 *
 *   dynamics.c  system
 *
 * with
 *
 *   cost.c incremental cost
 *
 *	first order descent (or 1.4) is the default
 *	second order is selected by defining
 *	#define SECOND_ORDER_DESCENT
 *
 *		** runs BACKWARD in time
 *		(in matlab, we must reverse the time trajectories ... :(  ))
 *		(and eat the '-' sign in -Pdot = A' P ... etc. )
 *		[latest matlab can probably go backwards, but stay the same]
 *
 * John Hauser
 * Dec 04
 * ...
 * Dec 11
 *
 *  states
 *	P - quadratic part of cost to go	n(n+1)/2
 *	r - linear part of cost to go		n
 *	q - *lagrange* multiplier for second order terms of projection operator
 *
 *  inputs (in time reversed order)
 *      xi(NS+NI)     alf(NS)   mu(NI) - trajectory for linearization
 *      K(NS*NI)      feedback gain for projection operator
 *      wt(NW)        exogenous system input w(t)
 *      wlt(NWL)      exogenous cost input w_c(t) (e.g., x_des, u_des)
 *
 *  outputs
 *      Ko(NS*NI)     new (optimal) feedback gain
 *                    Ko =   Ro^{-1} ( B^T P + S^T )
 *      vo(NI)        optimal feedforward input
 *                    vo = - Ro^{-1} ( B^T r + b )
 *  parameters
 *	none  ... ??? param(1)  [ Mass, cD0 ]
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
// NW  - number of exogenous inputs for dynamics
// NWL  - number of exogenous inputs for cost (& constraints ??)
// NPD  - number of parameters for dynamics
// NPC  - number of parameters for cost (& constraints ??)

#define NS2             (NS*NS)
#define NN12            ((NS*(NS+1))/2)

#define N_STATES        (NN12 + NS + NS)
#define N_INPUTS        (NS+NI + NS*NI + NW + NWL)
#define N_OUTPUTS       (NS*NI + NI)
#define N_PARAMS        (0)

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

static void optBKv(double *B_, double *Ko_, double *vo,
		  double *x, double *u)
{

#define B(i,j)	B_[ ((i)-1) + ((j)-1)*(NS) ]	// columnwise
#define Ko(i,j)	Ko_[ ((i)-1)*NS + ((j)-1) ]	// rowwise
#define P(i,j)	(( (j<=i) ? P_[i-1][j-1] : P_[j-1][i-1] ))	// low tri
#define S(i,j)	S_[ ((i)-1) + ((j)-1)*(NS) ]	// columnwise
#define R(i,j)  R_[ ((i)-1) + ((j)-1)*(NI) ]

#define q_fxu_x_u(i,j)  q_fxu_x_u_[((i)-1) + NS*((j)-1)]
#define q_fxu_u_u(i,j)  q_fxu_u_u_[((i)-1) + NI*((j)-1)]

#define lxu_x_u(i,j)    lxu_x_u_[((i)-1) + NS*((j)-1)]
#define lxu_u_u(i,j)    lxu_u_u_[((i)-1) + NI*((j)-1)]


  double *alf = u, *mu = alf+NS,
    *K_ = mu+NI, *wt = K_+NS*NI, *wlt = wt+NW;

  double *P__ = x, *r = P__+NN12, *q = r+NS;

  double *P_[NS];    // array of pointers for lower triangular P

  double Ri_[NI];    // for inverse of the *diagonal* R

  double  lxu_x_u_[NS*NI], lxu_u_u_[NI*NI];
  double *S_ = lxu_x_u_,  *R_ = lxu_u_u_;
  double b[NI];

  double q_fxu_x_u_[NS*NI], q_fxu_u_u_[NI*NI];

  int i, j, k;
  double ko;

  j = 0;	// set up lower triangular structure for P
  for (i=0; i<NS; i++) {
    P_[i] = P__ + j;
    j += i+1;
  }

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  // get derivs of the vector field and cost about  (alf(.),mu(.))
  // need B b S R

  // get b(4), S(16), R(32)
  cost(alf,mu,wlt,4+16+32, NULL, NULL, b,    NULL,S_,R_);

#ifdef SECOND_ORDER_DESCENT
  // get B(8), S(32), & R(64)
dynamics(alf,mu,wt, 8+32+64, NULL,NULL, NULL,B_, q, NULL,q_fxu_x_u_,q_fxu_u_u_);
#else
  // get B(8)
  dynamics(alf,mu,wt, 8, NULL,NULL, NULL,B_, NULL, NULL,NULL,NULL);
#endif

#ifdef SECOND_ORDER_DESCENT
  // add in vector field terms
  for (i=1; i<=NS; i++) {
    for (j=1; j<=NI; j++){
      S(i,j) += q_fxu_x_u(i,j);
    }
  }
#else
  // nothing to do
#endif // SECOND_ORDER_DESCENT

#ifdef SECOND_ORDER_DESCENT
  // add in vector field terms
  for (i=1; i<=NI; i++) {
    for (j=1; j<=NI; j++){
      R(i,j) += q_fxu_u_u(i,j);
    }
  }
#else
  // nothing to do
#endif

  // compute R^{-1} only for *diagonal* R ******
  //   also, we really need to check invertibility
  for (i=1; i<=NI; i++) {
    Ri_[i-1] = 1.0/R(i,i);
  }

  // Kopt = R^{-1} (B^T P + S^T)
  //     this calc requires R to be diagonal
  for (i=1; i<=NI; i++) {
    for (j=1; j<=NS; j++) {
      ko = 0.0;
      for (k=1; k<=NS; k++) {
	ko += B(k,i)*P(k,j);    // B^T P
      }
      Ko(i,j) = Ri_[i-1]*(ko + S(j,i));   // R^{-1} (B^T P + S^T)
    }
  }

  /* v_o =  - Ro^{-1} ( B^T r + b )  */

  for (i=1; i<=NI; i++) {
    vo[i-1] = - Ri_[i-1] * b[i-1];
    for (j=1; j<=NS; j++) {
      vo[i-1] -= Ri_[i-1]*B(j,i)*r[j-1];
    }
  }

#undef B
#undef Ko
#undef P
#undef S
#undef R

#undef q_fxu_x_u
#undef q_fxu_u_u

#undef lxu_x_u
#undef lxu_u_u
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
  real_T	*y	= ssGetOutputPortRealSignal(S,0);
  real_T	*x	= ssGetContStates(S);
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double     *param0 = mxGetPr(ssGetSFcnParam(S,0));

  int i;
  double u[N_INPUTS];
  double B_[NS*NI], *Ko_ = y, *vo = Ko_ + NS*NI;

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

  /* output new optimal gain and feedforward input */

  optBKv(B_, Ko_, vo, x, u);
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
  // double        *param0  = mxGetPr(ssGetSFcnParam(S,0));

  int i, j, k;
  double u[N_INPUTS];
  double *alf = u, *mu = alf+NS, *Kr_ = mu+NI,
    *wt = Kr_+NI*NS, *wlt = wt+NW;

  double a[NS], b[NI], Ko_b[NS], Kr_b[NS];

  double  *P__ = x,   *r = P__+NN12,   *q = r+NS;
  double *dP__ = dx, *dr = dP__+NN12, *dq = dr+NS;
  double dp;
  double *P_[NS], *dP_[NS];	// array of pointers for lower triangular P, dP

  double A_[NS2], A_BKo_[NS2], A_BKr_[NS2];
  double B_[NS*NI], Ko_[NI*NS];
  double Q_[NS2], R_[NI*NI];

  double q_fxu_x_x_[NS2], q_fxu_u_u_[NI*NI];

  double vo[NI];

  /* index A & P with 1:NS  (NOT 0:NS-1!) */
#define A(i,j)        A_[((i)-1) + NS*((j)-1)]
#define A_BKo(i,j)    A_BKo_[((i)-1) + NS*((j)-1)]
#define A_BKr(i,j)    A_BKr_[((i)-1) + NS*((j)-1)]
#define Q(i,j)        Q_[((i)-1) + NS*((j)-1)]
#define R(i,j)        R_[((i)-1) + NI*((j)-1)]
#define P(i,j)        (( (j<=i) ? P_[i-1][j-1]  : P_[j-1][i-1] ))  // low tri
#define dP(i,j)       ( dP_[i-1][j-1] )           // *** i <= j *** REQUIRED!
#define B(i,j)        B_[ ((i)-1) + ((j)-1)*(NS) ]      // columnwise
#define Kr(i,j)       Kr_[ ((i)-1)*NS + ((j)-1) ]     // rowwise
#define Ko(i,j)       Ko_[ ((i)-1)*NS + ((j)-1) ]     // rowwise

#define q_fxu_x_x(i,j)  q_fxu_x_x_[((i)-1) + NS*((j)-1)]
#define q_fxu_u_u(i,j)  q_fxu_u_u_[((i)-1) + NI*((j)-1)]

  // set up lower triangular structure for P and dP
  j = 0;
  for (i=0; i<NS; i++) {
    P_[i]  = P__ + j;
    dP_[i] = dP__ + j;
    j += i+1;
  }

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  // get 1st & 2nd derivs of vector field and cost about  (alf(.),mu(.))
  // need a b A Q R

  // get a(2), b(4), Q(8), R(32)
  cost(alf,mu,wlt,2+4+32+64, NULL, a, b,    Q_,NULL,R_);

#ifdef SECOND_ORDER_DESCENT
  // get A(4), Q(16), & R(64)
dynamics(alf,mu,wt, 4+16+64, NULL,NULL, A_,NULL, q, q_fxu_x_x_,NULL,q_fxu_u_u_);
#else
  // get A(4)
  dynamics(alf,mu,wt, 4,      NULL,NULL, A_,NULL, NULL, NULL,NULL,NULL);
#endif

  optBKv(B_, Ko_, vo, x, u);

#ifdef SECOND_ORDER_DESCENT
  for (i=1; i<=NS; i++){
    for (j=1; j<=NS; j++){
      Q(i,j) += q_fxu_x_x(i,j);
    }
  }
  for (i=1; i<=NI; i++){
    for (j=1; j<=NI; j++){
      R(i,j) += q_fxu_u_u(i,j);
    }
  }
#else
  // nothing
#endif //SECOND_ORDER_DESCENT

  /* only *lower triangle* of dP */
  for (i=1; i<=NS; i++) {
    for (j=1; j<=i; j++) {
      dp = Q(i,j);
      for (k=1; k<=NS; k++) {
        dp += A(k,i)*P(k,j) + P(i,k)*A(k,j);    // A^T P + P A
      }
      for (k=1; k<=NI; k++) {
        // assume R has only terms on diagonal here
        dp -= Ko(k,i)*R(k,k)*Ko(k,j);   // -P B R^-1 B^T P = -Ko^T R Ko
      }
      dP(i,j) /* = dP(j,i) */  = dp;
    }
  }

  /* make copies of A(.) to use with K(.) and K_opt(.) */
  for (i=0; i<NS2; i++)
    A_BKr_[i] = A_BKo_[i] = A_[i];

  /* feedbacks, optimal and projection */
  for (i=1; i<=NS; i++) {
    for (j=1; j<=NS; j++) {
      for (k=1; k<=NI; k++) {
        A_BKo(i,j) -= B(i,k)*Ko(k,j);
        A_BKr(i,j) -= B(i,k)*Kr(k,j);
      }
    }
  }

  // K^T b
  for (i=1; i<=NS; i++) {
    Kr_b[i-1] = Ko_b[i-1] = 0.0;
    for (j=1; j<=NI; j++) {
      Ko_b[i-1] += Ko(j,i)*b[j-1];
      Kr_b[i-1] += Kr(j,i)*b[j-1];
    }
  }

  /* "r" uses the optimal feedback gain Ko */
  /* "q" depends on the projection operator feedback gain K (or Kr) */

  for (i=1; i<= NS; i++) {
    dr[i-1] = 0.0 + a[i-1] - Ko_b[i-1];
    dq[i-1] = 0.0 + a[i-1] - Kr_b[i-1];
    for (j=1; j<=NS; j++) {
      dr[i-1] += A_BKo(j,i)*r[j-1];
      dq[i-1] += A_BKr(j,i)*q[j-1];
    }
    // printf("dq[%d]: %g\n", i, dq[i-1]);
  }

#undef Kr
#undef Ko
#undef B
#undef dP
#undef P
#undef R
#undef Q
#undef A_BKr
#undef A_BKo
#undef A

#undef q_fxu_x_x
#undef q_fxu_u_u
}

static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

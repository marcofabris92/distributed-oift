/*
 * dynamics.c
 *
 *
 *    for Formation Control and Path Following
 *
 *      dx = Ax + Bu
 *
 * dynamics(x, u, wt, ders, dx,    y,    fxu_x_, fxu_u_, q, q_fxu_x_x_, q_fxu_x_u_, q_fxu_u_u_);
 * ders:                    dx(1), y(2), A(4),   B(8),      Q(16),      S(32),      R(64)
 *       q  is the 'projection operator stabilized adjoint q'
 *
 * state             (x)
 *
 *   p1, ... , pn, dp1, ... , dpn
 *
 * input             (u)
 *
 *   u1, ..., un
 *
 * exogenous input   (wt)
 *
 *   x_des (pB, dpB)
 *
 * output            (y)
 *
 *   none
 *
 *
 * Marco Fabris, Boulder, 06/29/2018
 * checked out
 *
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wct,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),  Q(8),Qsafe(64), S(16), R(32)
 */

#include <math.h>
#include "sys_sizes.h"

/*
 * dynamics(x, u, wt, ders, dx,    y,    fxu_x_, fxu_u_, q, q_fxu_x_x_, q_fxu_x_u_, q_fxu_u_u_);
 * ders:                    dx(1), y(2), A(4),   B(8),      Q(16),      S(32),      R(64)
 */

void dynamics(
  double *x, double *u, double *wt,
  int ders,
  double *dx, double *y,
  double *fxu_x_, double *fxu_u_,
  double *q,
  double *q_fxu_x_x_, double *q_fxu_x_u_, double *q_fxu_u_u_
  )
{

#define NX  (NS)
#define NU  (NI)

#define A(i,j)   fxu_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define B(i,j)   fxu_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise

#define q_fxu_x_x(i,j)  q_fxu_x_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define q_fxu_x_u(i,j)  q_fxu_x_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define q_fxu_u_u(i,j)  q_fxu_u_u_[ ((i)-1) + ((j)-1)*(NU) ]  // columnwise

  int do_dx, do_A, do_B, do_Q, do_S, do_R, do_y;

  int i, j, k;

  // determine what jobs to do
  /*
   *    ders - which ders? use binary code to specify
   *    1 - fxu
   *    2 - y
   *    4 - fxu_x = A
   *    8 - fxu_u = B
   *   16 - q_fxu_x_x = Q
   *   32 - q_fxu_x_u = S
   *   64 - q_fxu_u_u = R
   */
  do_dx =  ders & 1       ;
  do_y  = (ders & 2)  >> 1;
  do_A  = (ders & 4)  >> 2;
  do_B  = (ders & 8)  >> 3;
  do_Q  = (ders & 16) >> 4;
  do_S  = (ders & 32) >> 5;
  do_R  = (ders & 64) >> 6;


  // UPDATING STATE SPACE SYSTEM

  if (do_A || do_Q || do_S) {
    // nothing here
  }

  if (do_A) {
    // zero out A
    for (i=0; i<NX*NX; i++) {
      fxu_x_[i] = 0.0;
    }
  }

  if (do_B) {
    // zero out B
    for (i=0; i<NX*NU; i++) {
      fxu_u_[i] = 0.0;
    }
  }

  if (do_Q) {
    // zero out Q = q_fxu_x_x
    for (i=0; i<NX*NX; i++) {
      q_fxu_x_x_[i] = 0.0;
    }
  }

  if (do_S) {
    // zero out S = q_fxu_x_u
    for (i=0; i<NX*NU; i++) {
      q_fxu_x_u_[i] = 0.0;
    }
  }

  if (do_R) {
    // zero out R = q_fxu_u_u
    for (i=0; i<NU*NU; i++) {
      q_fxu_u_u_[i] = 0.0;
    }
  }

  if (do_dx) {

    for(i = 0; i < NU; i++){
    	dx[i] = x[i+NU];
    	dx[i+NU] = u[i];
    }

  }

  if (do_y) {

    // no term in here

  }

  if (do_A) {

  	for(i = 1; i <= NU; i++)
  		A(i,i+NU) = 1.0;

  }

  if (do_B) {

    for(i = 1; i <= NU; i++)
  		B(i+NU,i) = 1.0;

  }


  if (do_Q) {

  	// no term in here

}

  if (do_S) {

	// no term in here

  }

  if (do_R) {

    // no term in here

  }

#undef A
#undef B

#undef q_fxu_x_x
#undef q_fxu_x_u
#undef q_fxu_u_u

#undef NX
#undef NU
}



/*
 *  dynamics_unbdd(x,u,wt, dx);
 *
 *  trap to freeze a system that is becoming unbounded
 */
 void dynamics_unbdd(
  double *x, double *u, double *wt,
  double *dx
  )
{
  /*
#define NX (NS)
#define NU (NI)

  int i;

  // if pendulum angle exceeds 30 radians, just quit! 		// do we need to change this?
  if ( x[0]*x[0] > 900.0 ) {
    for (i=0; i<NX; i++) {
      dx[i] = 0.0;
    }
  }

#undef NX
#undef NU
*/
}

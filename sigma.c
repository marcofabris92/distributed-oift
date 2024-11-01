#include <math.h>

// adapted from  sig.m  :
//
// function [sg, dsig, ddsig] = sig(s,d,k)
//
//      s can be a vector of nonnegative values
//        that should be appropriately restricted depending
//        on which function is desired and what derivs are needed
//
//      d is a positive scalar
//
//      k = 0 or 1
//
//      sig0(s) = c0*(1-s/(d^2))^3																	     s < d^2
//
//      sig1(s) = c1*( sqrt( s )           - d )^3 / d^3     						 s > d^2
//
//    which we hope can be glued together to give a C^2 function on s >= 0 .
//
//    note that both of these functions are well defined for  0 <= s <= 2*d^2
//    and smooth (even analytic) on  0 < s < 2*d^2
//
//    thus we can compare the two functions (+ derivs) on that domain
//
//
//  JH jul18 boulder, MF oct06 camponogara

double sigma(double s, double d, double *sig_s, double *sig_s_s) {
	double d2;
  double scale0, scale1;
  double sig, dsig, ddsig;

  if (s < 0.0) { s = -s; }

  if (d < 0.0) { d = -d; }

  d2 = d*d;
  scale0 = 250.0;  // 100 when dd = 1 for each inter-agent distance, 1 when shaping is required --> kr
  scale1 = 100.0/(d*d2); // 1 when dd = 1 for each inter-agent distance, 10 when shaping is required --> ka

  if ( s >= d2 ) {
  	double r = sqrt(s);
    double r_d = r-d;
    double r_d2 = r_d*r_d;

    // sig = ( sqrt( s ) - d ).^3;    						% sig1
    sig = r_d2*r_d;
    sig *= scale1;

    // dsig = 1.5*((sqrt(s) -d ).^2)./sqrt(s);    % dsig1
    dsig = 1.5*r_d2/r;
    dsig *= scale1;

    // ddsig = 0.75*(s-d^2)/(r^3);    						% ddsig1
    ddsig = 0.75*(s-d2)/(s*r);
    ddsig *= scale1;

  } else {

  	double base = 1.0-s/d2;
  	double base2 = base*base;

    // sg = (1-s/d2)^3;  																					% sig0
    sig = base2*base;
    sig *= scale0;

    // dsig = (-3/d2)*(1-s/d2)^2;    															% dsig0
    dsig = -3.0/d2*base2;
    dsig *= scale0;

    // ddsig = 6/(d2*d2)*(1-s/d2) 																% ddsig0
    ddsig  = 6.0/(d2*d2)*base;
    ddsig *= scale0;
  }

  *sig_s   = dsig;
  *sig_s_s = ddsig;
  return sig;
}

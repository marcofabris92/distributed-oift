/*
 *  cost_params.h
 *
 *  for Formation Control and Path Following
 */

#include "sys_sizes.h"


// QB    | 2*DIMx2*DIM
// RD    | NIxNI
// SD = 0| NSxNI
// P1    | NSxNS

// diagonal QB
#define q_p (10.0/nAg) // 10 for complex tracking // 10 for normal tracking
#define q_v (1.0/nAg)

// diagonal RD
#define r_a (1.0)

// formation weight for FF term
#define kF (2.0) // 1 0.4

// velocity alignment weight for AA term
#define kA (0.25) // 1


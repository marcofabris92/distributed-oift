/*
 * sys_sizes.h
 *
 * specify the system sizes (NS, NI, etc.)
 *
 * this should be replicated in
 *
 *   sys_sizes_h.m
 *
 */

// Formation Control and Path Following

#define dd    (5.0)				 // desired distance to keep in formation
#define DIM   (3)		       // dimension of the space in which agents live
#define nAg   (8)		       // number of agents
#define NS    (2*DIM*nAg)  // number of system states
#define NI    (DIM*nAg)    // number of system inputs
#define NO    (0)          // number of system outputs
#define NW    (0)          // number of exogenous inputs for dynamics
#define NWL   (2*DIM)      // number of exogenous inputs for cost l(x,u,w_l(t))
#define NWC   (0)          // number of constraints (0 if no constraints)
#define DT		(0.01)			 // dt used in the main script for the time discretization

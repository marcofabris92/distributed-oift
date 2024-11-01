# distributed-oift
This code replicates the numerical simulations reported in the research paper entitled "Optimal time-invariant distributed formation tracking for second-order multi-agent systems"

*********************************************************************************************************************************************************************

Authors: Marco Fabris*, Giulio Fattore, Angelo Cenedese

All the authors are with the University of Padua, Italy
 
\* M. Fabris is the algorithm and software developer. E-mail: marco.fabris.1@unipd.it
  
Special thanks to John Hauser** for his assistance while using PRONTO.

\** John Hauser is with the University of Colorado-Boulder, USA

Published on the European Journal of Control. 
Paper available at https://doi.org/10.1016/j.ejcon.2024.100985
 
Publication history:
- Received 23 June 2023
- Revised 31 December 2023
- Accepted 27 March 2024
- Available online 4 April 2024
- Version of Record 8 April 2024


Abstract:
This paper addresses the optimal time-invariant formation tracking 
problem with the aim of providing a distributed solution for multi-agent 
systems with second-order integrator dynamics. In the literature, most of
the results related to multi-agent formation tracking do not consider 
energy issues while investigating distributed feedback control laws. In 
order to account for this crucial design aspect, we contribute by 
formalizing and proposing a solution to an optimization problem that 
encapsulates trajectory tracking, distance-based formation control and 
input energy minimization, through a specific and key choice of potential
functions in the optimization cost. To this end, we show how to compute 
the inverse dynamics in a centralized fashion by means of the 
Projector-Operator-based Newton’s method for Trajectory Optimization 
(PRONTO) and, more importantly, we exploit such an offline solution as a 
general reference to devise a stabilizing online distributed control law. 
Finally, numerical examples involving a cubic formation following a 
chicane-like path in the 3D space are provided to validate the proposed 
control strategies.
*********************************************************************************************************************************************************************

Instructions:

0. install a GCC compiler on your PC (see also https://it.mathworks.com/help/matlab/call-mex-files-1.html?s_tid=CRUX_lftnav) and configure it by running mex_newt_files.m in MatLab (make sure to indicate the correct path in the command setenv)
1. set the kind of simulation you wish to replicate by assigning the variable simul_choice in MAIN.m (or create a new simulation)
2. run MAIN.m and wait for its termination


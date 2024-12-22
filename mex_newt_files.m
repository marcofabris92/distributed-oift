setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
% mex -setup

% make dynamics and cost functions accessible as matlab commands
mex dynamics_m.c
mex cost_m.c
mex sys_sizes_m.c
% mex QR_params_m.c

% regulator
% mex LQR_Kr.c

% projection operator & cost
mex nonl_K.c

% descent direction calcs
mex Prq1_Kv.c
mex Prq2_Kv.c
mex linK_zv.c


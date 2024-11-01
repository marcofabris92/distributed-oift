function [xxdot] = distr_dyn(t,xx,par)


% state
% x = [position   velocity]

global flag33 flag50 flag90

dt = par.dt;
nAg = par.nAg;
M = par.DIM;
QBp = par.QBp;
QBdp = par.QBdp;
R = par.R;
kF = par.kF;
kA = par.kA;
qA = par.qA;
Xdes = par.Xdes;
dijs = par.dijs;
tl = par.tl;
kr = par.kr;
ka = par.ka;
G = par.G;
T0 = par.T0;
dist_coeff = par.dist_coeff;
max_dist = par.max_dist;
gains = par.gains;
gammaOBS_P = par.gammaOBS_P;
gammaOBS_D = par.gammaOBS_D;
tt = 1+floor(t/dt);
NI = nAg*M;
NS = NI*2;
xdot = zeros(NS,1);


if t/tl > 0.33 && ~flag33
    fprintf('> 33%%\n')
    flag33 = 1;
end
if t/tl > 0.50 && ~flag50
    fprintf('> 50%%\n')
    flag50 = 1;
end
if t/tl > 0.90 && ~flag90
    fprintf('> 90%%\n')
    flag90 = 1;
end


%% integrating the dynamics

pBdes = Xdes(tt,1:M)';
dpBdes = Xdes(tt,M+1:M+M)';

p_hat = zeros(NI,nAg);
dp_hat = zeros(NI,nAg);


% velocity integrator ---> pdot = v
xdot(1:nAg*M) = xx(nAg*M+1:2*nAg*M);

% true centroid (just in case we need to compare performances)
pB = zeros(M,1);
dpB = zeros(M,1);
for i = 1:nAg
    intval_i = (i-1)*M+1:(i-1)*M+M;
    pB = pB + xx(intval_i)/nAg;
    dpB = dpB + xx(intval_i+M*nAg)/nAg;
end

% input to be assigned to the velocity ---> vdot = u
I_M = eye(M);
kP_tr = gains(1);
kD_tr = gains(2);
kP_fo = gains(3);
kD_fo = gains(4);
kP_al = gains(5);
NI2 = NI*nAg;
for i = 1:nAg
    Ni = neighbors(G,i);
    deg_i = length(Ni);
    formation_i = zeros(M,1);
    alignment_i = zeros(M,1);
    dformation_i = zeros(M,1);
    intval_i = (i-1)*M+1:(i-1)*M+M;
    p_i = xx(intval_i);
    dp_i = xx((nAg+i-1)*M+1:(nAg+i-1)*M+M);
    p_hat(:,i) = xx(NS+(i-1)*NI+1:NS+(i-1)*NI+NI);
    dp_hat(:,i) = xx(NI2+NS+(i-1)*NI+1:NI2+NS+(i-1)*NI+NI);
    for jj = 1:deg_i
        j = Ni(jj);
        dij = dijs(i,j);
        p_j = xx((j-1)*M+1:(j-1)*M+M);
        dp_j = xx((nAg+j-1)*M+1:(nAg+j-1)*M+M);
        eij = p_i-p_j;
        deij = dp_i-dp_j;
        sij = eij'*eij;
        sig1_ij = sigma(sij,dij,1,kr,ka);
        sig2_ij = sigma(sij,dij,2,kr,ka);
        sig1_h = sig1_ij > 0;
        if 0
            sig1_h = 1;
        end
        newt_ij = 1; 
        formation_i = formation_i + kF*newt_ij*sig1_ij*eij;
        alignment_i = alignment_i + kA*qA(i,j)*deij;  
        dformation_i = dformation_i + kF*(2*sig2_ij*(eij*eij')+...
            sig1_h*sig1_ij*I_M)*deij;
    end
    
    % acceleration integrator ---> pddot = u
    pBhat = zeros(M,1);
    dpBhat = zeros(M,1);
    for k = 1:nAg
        pBhat = pBhat + p_hat((k-1)*M+1:(k-1)*M+M,i)/nAg;
        dpBhat = dpBhat + dp_hat((k-1)*M+1:(k-1)*M+M,i)/nAg;
    end
    
    % pretending to use perfect information if we wish
    if 0
        pBhat = pB;
        dpBhat = dpB;
    end
    
    %E = [norm(pBhat-pBdes)^2 norm(dpBhat-dpBdes)^2];
    %E
    
    %% control law
    xdot(nAg*M+intval_i) = -(R(intval_i,intval_i)^-1)*...
        (kP_tr*QBp*(pBhat-pBdes) +...
         kD_tr*QBdp*(dpBhat-dpBdes) +...
         kP_fo*formation_i +...
         kD_fo*dformation_i +...
         kP_al*alignment_i);
    
     %[norm(mean(p_hat(:,i))-pBdes) norm(mean(dp_hat(:,i))-dpBdes)]
     %pause

end



%% checking saturation
u_sat = 50;
flag_t = 0;
for j = 1:nAg*M
    u_j = xdot(nAg*M+j);
    if abs(u_j) > u_sat
        xdot(nAg*M+j) = sign(u_j)*u_sat;
        if ~flag_t
            fprintf(strcat('Saturation for t =',num2str(t),'\n'))
            flag_t = 1;
        end
    end
end
% in this example the artificial saturation works from t = 0
% to t = 0.35726, since the initial conditions for u is inappropriate

    

%% observer for the centroid
up_obs = zeros(2*NI2,1);
for i = 1:nAg
    Ni = neighbors(G,i);
    intval_i = (i-1)*NI+1:(i-1)*NI+NI;
    idp = zeros(NI,1);
    iu = zeros(NI,1);
    intval_ii = (i-1)*M+1:(i-1)*M+M;
    idp(intval_ii) = dp_hat(intval_ii,i); %xdot(intval_ii);
    iu(intval_ii) = xdot(NI+intval_ii);
    for j_ = 1:length(Ni)
        j = Ni(j_);
        intval_jj = (j-1)*M+1:(j-1)*M+M;
        up_obs(intval_i) = up_obs(intval_i) +...
            p_hat(:,i)-p_hat(:,j);
        up_obs(NI2+intval_i) = up_obs(NI2+intval_i) +...
            dp_hat(:,i)-dp_hat(:,j);
        idp(intval_jj) = dp_hat(intval_jj,i); %xdot(intval_jj);
        iu(intval_jj) = xdot(NI+intval_jj);
    end
    notNi = setdiff(1:nAg,[i; Ni]);
    for j_ = 1:length(notNi)
        j = Ni(j_);
        intval_jj = (j-1)*M+1:(j-1)*M+M;
        idp(intval_jj) = dp_hat(intval_jj,i);
    end

    projp = zeros(NI,1);
    projdp = zeros(NI,1);
    projp(intval_ii) = p_hat(intval_ii,i)-xx(intval_ii);
    projdp(intval_ii) = dp_hat(intval_ii,i)-xx(NI+intval_ii);
    up_obs(intval_i) = -gammaOBS_P*(up_obs(intval_i) + projp) + idp;
    up_obs(NI2+intval_i) = -gammaOBS_D*(up_obs(NI2+intval_i) + projdp) +...
        iu;  
end



%% final update
% pdot = xdot(1:NI);
% dpdot = xdot(NI+1:NS);
xxdot = [xdot; up_obs];




    
end


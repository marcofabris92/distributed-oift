close all
clear all, clc

global dt nAg M kr ka QBp QBdp R kF kA Xdes dijs G flag33 flag50 flag90

dt = 0.01;
tl = 200;
T0 = (0:dt:tl)';
T = length(T0);

nAg = 4;
I_nAg = eye(nAg);
M = 2;
I_M = eye(M);

QBp = 10^1;
QBdp = 10^1;
R = 10^0;
kF = 10^0;
kA = 5*10^-1;
kr = 1;
ka = 50;

% initial conditions
vel = 1;
x0p = [-2 1  -3 -1  2 -2  0 0];   
x0p = x0p/(2*norm(x0p));
x0p = rand(1,8)*10;
x0v = 5*[0 -vel  0 -vel  0 -vel  0 -vel];

% desired trajectory for the centroid
xd0p = zeros(1,M);
xd0v = 1*[vel 0];
Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
Xdesv = repmat(xd0v,length(T0),1);
Xdes = [ Xdesp Xdesv ];

% framework definition
dist_coeff = 1.05;
dd = 5;
% dijs = [0 dd dd;
%         dd 0 dd;
%         dd dd 0];
dd1_2 = sqrt(2)*dd;
max_dijs = dd1_2;
INF = eps+dist_coeff*dd1_2;
dijs = [0 dd INF dd;
        dd 0 dd dd1_2;
        INF dd 0 dd;
        dd dd1_2 dd 0];

A = zeros(nAg);
for i = 1:nAg
    for j = 1:i-1
        if dijs(i,j) > 0
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end
G = graph(A);
degs = sum(A);
D = diag(degs);
L = D-A;
L_ = D^(-1/2)*L*D^(-1/2);
eigsL_ = sort(eig(L_));
lambda1 = eigsL_(2);
lambdan_1 = eigsL_(end);
muL_ = (lambda1+lambdan_1)/2;
eta_star = 1-1/muL_;
if eta_star < 0
    eta_star = 0;
end
% if topology remains constant:
F = eye(nAg)*eta_star+(D^-1*A)*(1-eta_star);


%% -------------------------------------------------------
% parameters
params.dt = dt;
params.nAg = nAg;
params.M = M;
params.QBp = QBp;
params.QBdp = QBdp;
params.R = R;
params.kF = kF;
params.kA = kA;
params.Xdes = Xdes;
params.dijs = dijs;
params.G = G;
params.tl = tl;
params.F = kron(F,I_M);
params.dist_coeff = dist_coeff;
params.max_dist = max_dijs; % max(dijs(:));

% flags to show completion
flag33 = 0;
flag50 = 0;
flag90 = 0;

% dynamics integration
[t,x_history] = ode45(@(t,x)distr_dyn(t,x,params),T0',[x0p x0v]');


%% display final results
final_state = x_history(end,:)';
disp(final_state)
p1 = x_history(end,1:2)'
p2 = x_history(end,3:4)'
p3 = x_history(end,5:6)'
p4 = x_history(end,7:8)'
e12 = p1-p2;
e13 = p1-p3;
e23 = p2-p3;
e14 = p1-p4;
e24 = p2-p4;
e34 = p3-p4;
pB = (p1+p2+p3+p4)/4
Ne12 = norm(e12)
Ne13 = norm(e13)
Ne23 = norm(e23)
Ne14 = norm(e14)
Ne24 = norm(e24)
Ne34 = norm(e34)


%% fase diagram for positions



figure
grid on
hold on
plot(Xdes(:,1),Xdes(:,2),'k')
for i = 1:nAg
    plot(x_history(:,(i-1)*M+1),x_history(:,(i-1)*M+2),'r')
end
for i = 1:nAg
    plot(x_history(end,(i-1)*M+1),x_history(end,(i-1)*M+2),'b*','linewidth',3)
end
% inserire qui un confronto tra grafici... prendere i dati altrui potrebbe
% essere una buona idea per poi plottare la traiettoria distribuita su
% quella centralizzata -> fare un esempio con triangolo equilatero
% - traiettoria
% - funzionale di costo
% - input speso
% - ripartizione tra formation vs tracking
% * consensus solo per il distribuito
% nota: mettere il peso finale per la final boundary condition nel caso
% centralizzato!! (peso finale appropriato, no un punto a caso per la
% traiettoria desiderata)
axis equal

%% evolutions of the distance errors (sigma_dij-s)

figure
grid on
hold on
dist_errors = zeros(T,nAg*(nAg-1)/2);

k = 0;
for i = 2:nAg
    for j = 1:i-1
        dij = dijs(i,j);
        if dij >= 0
            for tt = 1:T
                p_i_t = x_history(tt,(i-1)*M+1:(i-1)*M+M);
                p_j_t = x_history(tt,(j-1)*M+1:(j-1)*M+M);
                sij_t = norm(p_i_t-p_j_t)^2;
                dist_errors(tt,k+1) = sigma(sij_t,dij,0);
            end
        else
            dist_errors(:,k+1) = -ones(T,1);
        end
        k = k+1;
    end
end

k = 0;
for i = 2:nAg
    for j = 1:i-1
        tratto = '-';
        if dijs(i,j) < 0
            tratto = '--';
        end
        plot(T0,dist_errors(:,k+1),tratto,'linewidth',1.5)
        k = k+1;
    end
end

title('zero-order consensus')
xlabel('$t$','interpreter','latex')
ylabel('$\sigma(s_{ij})\qquad$','interpreter','latex')
set(gca,'fontsize',25)
set(get(gca,'ylabel'),'rotation',0)


%% evolutions of the first derivatives of the distance errors (sigma_dij-s)

figure
grid on
hold on
ddist_errors = zeros(T,nAg*(nAg-1)/2);

k = 0;
for i = 2:nAg
    for j = 1:i-1
        dij = dijs(i,j);
        if dij >= 0
            for tt = 1:T
                p_i_t = x_history(tt,(i-1)*M+1:(i-1)*M+M);
                p_j_t = x_history(tt,(j-1)*M+1:(j-1)*M+M);
                sij_t = norm(p_i_t-p_j_t)^2;
                ddist_errors(tt,k+1) = sigma(sij_t,dij,1);
            end
        else
            ddist_errors(:,k+1) = zeros(T,1);
        end
        k = k+1;
    end
end

k = 0;
for i = 2:nAg
    for j = 1:i-1
        tratto = '-';
        if dijs(i,j) < 0
            tratto = '--';
        end
        plot(T0,ddist_errors(:,k+1),tratto,'linewidth',1.5)
        k = k+1;
    end
end

title('first-order consensus')
xlabel('$t$','interpreter','latex')
ylabel('$\frac{\partial\sigma}{\partial s_{ij}}\qquad$','interpreter','latex')
set(gca,'fontsize',25)
set(get(gca,'ylabel'),'rotation',0)

%% evolutions of the second derivatives of the distance errors (sigma_dij-s)

figure
grid on
hold on
dddist_errors = zeros(T,nAg*(nAg-1)/2);

k = 0;
for i = 2:nAg
    for j = 1:i-1
        dij = dijs(i,j);
        if dij >= 0
            for tt = 1:T
                p_i_t = x_history(tt,(i-1)*M+1:(i-1)*M+M);
                p_j_t = x_history(tt,(j-1)*M+1:(j-1)*M+M);
                sij_t = norm(p_i_t-p_j_t)^2;
                dddist_errors(tt,k+1) = sigma(sij_t,dij,2);
            end
        else
            dddist_errors(:,k+1) = zeros(T,1);
        end
        k = k+1;
    end
end

k = 0;
for i = 2:nAg
    for j = 1:i-1
        tratto = '-';
        if dijs(i,j) < 0
            tratto = '--';
        end
        plot(T0,dddist_errors(:,k+1),tratto,'linewidth',1.5)
        k = k+1;
    end
end

title('second-order consensus')
xlabel('$t$','interpreter','latex')
ylabel('$\frac{\partial^2\sigma}{\partial s_{ij}^2}\qquad$','interpreter','latex')
set(gca,'fontsize',25)
set(get(gca,'ylabel'),'rotation',0)


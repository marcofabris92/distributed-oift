% "Optimal time-invariant distributed formation tracking 
% for second-order multi-agent systems"
% 
%
% Authors:
% Marco Fabris*, Giulio Fattore, Angelo Cenedese
% All the authors are with the University of Padua, Italy
% 
% * M. Fabris is the algorithm and software development.
% e-mail: marco.fabris.1@unipd.it
% 
% Special thanks to John Hauser** for his assistance while using PRONTO.
% ** John Hauser is with the University of Colorado-Boulder, USA
%
% Published on the European Journal of Control
% 
% Publication history:
% - Received 23 June 2023
% - Revised 31 December 2023
% - Accepted 27 March 2024
% - Available online 4 April 2024
% - Version of Record 8 April 2024
%
%
% Abstract:
% This paper addresses the optimal time-invariant formation tracking 
% problem with the aim of providing a distributed solution for multi-agent 
% systems with second-order integrator dynamics. In the literature, most of
% the results related to multi-agent formation tracking do not consider 
% energy issues while investigating distributed feedback control laws. In 
% order to account for this crucial design aspect, we contribute by 
% formalizing and proposing a solution to an optimization problem that 
% encapsulates trajectory tracking, distance-based formation control and 
% input energy minimization, through a specific and key choice of potential
% functions in the optimization cost. To this end, we show how to compute 
% the inverse dynamics in a centralized fashion by means of the 
% Projector-Operator-based Newton’s method for Trajectory Optimization 
% (PRONTO) and, more importantly, we exploit such an offline solution as a 
% general reference to devise a stabilizing online distributed control law. 
% Finally, numerical examples involving a cubic formation following a 
% chicane-like path in the 3D space are provided to validate the proposed 
% control strategies.





%global NS NI QB R nAg DIM T0 dijs kF kA kr ka ftsz lw dt QBp QBdp G...
global    flag33 flag50 flag90

DIM = 3; % = M-dimensional space
kTgif = 100;

% get system sizes: NS, NI, NO, etc.
[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;
nAg = NI/DIM; % = n agents
dt = 0.01;
vel = 1;

dd = 5;
dijs = [];
A = zeros(nAg);

max_dijs = 5;
T = 0;
x11 = zeros(1,NS);
Rot = @(theta) ([cos(theta) -sin(theta); sin(theta) cos(theta)]);
Rzx = 0.5*[sqrt(2) 1 1; -sqrt(2) 1 1; 0 -sqrt(2) sqrt(2)];
Rzy = 0.5*[1 -sqrt(2) -1; 1 sqrt(2) -1; sqrt(2) 0 sqrt(2)];


%% Simulations settings
simul_choice = 1;
switch DIM
    case 2
        switch simul_choice
            
            case 1 % pentagon without few edges
                pdi = 2*sin(54*pi/180);
                dijs = [0      dd     dd*pdi dd*pdi dd;
                        dd     0      dd     -1     dd*pdi;
                        dd*pdi dd     0      dd     -1;
                        dd*pdi -1     dd     0      dd;
                        dd     dd*pdi -1     dd     0      ];
                QBp = 10^1/nAg;
                QBdp = 10^0/nAg;
                ptilde_p = 10/nAg; % terminal QBp 
                ptilde_v = 1/nAg; % terminal QBdp
                QB = [QBp*eye(DIM) zeros(DIM);
                      zeros(DIM)   QBdp*eye(DIM)];
                R = 1*eye(nAg*DIM);
                kF = 2;
                kA = 0.25;
                kr = 250;
                ka = 100;
                T = 20;
                T0 = (0:dt:T)';
                x0p = [-2 1  -3 -1  2 -2  5 1  -1 -6];   % randomizer 1*(2*(rand(1,DIM*nAg)-0.5));
                x0p = x0p/(2*norm(x0p));
                x0v = 5*[0 -vel  0 -vel  0 -vel 0 -vel 0 -vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel 0];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                T0e = T0(end);
                x11p = [];
                th = 72*pi/180;
                for i = 1:nAg
                    x11p = [x11p vel*T0e+dd*cos((i-1)*th) dd*sin((i-1)*th)];
                end
                x11v = [vel 0 vel 0 vel 0 vel 0 vel 0]; % zeros(1,10);
                x11 = [ x11p x11v ];
                kT = 5;
                
            case 2 % full pentagon 
                pdi = 2*sin(54*pi/180);
                dijs = [0      dd     dd*pdi dd*pdi dd;
                        dd     0      dd     dd*pdi dd*pdi;
                        dd*pdi dd     0      dd     dd*pdi;
                        dd*pdi dd*pdi dd     0      dd;
                        dd     dd*pdi dd*pdi dd     0      ];
                QBp = 10^1/nAg;
                QBdp = 10^0/nAg;
                ptilde_p = 10/nAg; % terminal QBp 
                ptilde_v = 1/nAg; % terminal QBdp
                QB = [QBp*eye(DIM) zeros(DIM);
                      zeros(DIM)   QBdp*eye(DIM)];
                R = 1*eye(nAg*DIM);
                kF = 2;
                kA = 0.25;
                kr = 250;
                ka = 100;
                T = 20;
                T0 = (0:dt:T)';
                x0p = [-2 1  -3 -1  2 -2  5 1  -1 -6];   % randomizer 1*(2*(rand(1,DIM*nAg)-0.5));
                x0p = x0p/(2*norm(x0p));
                x0v = 5*[0 -vel  0 -vel  0 -vel 0 -vel 0 -vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel 0];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                T0e = T0(end);
                x11p = [];
                th = 72*pi/180;
                for i = 1:nAg
                    x11p = [x11p vel*T0e+dd*cos((i-1)*th) dd*sin((i-1)*th)];
                end
                x11v = [vel 0 vel 0 vel 0 vel 0 vel 0]; % zeros(1,10);
                x11 = [ x11p x11v ];
                kT = 5;
        end
        
  %% 3D simulations 
    case 3
        switch simul_choice
            
            case 1
                % cube without few edges
                sq2 = sqrt(2);
                sq3 = sqrt(3);
                dijs = dd*[0 1 -sq2 1 1 sq2 -sq3 sq2;
                           1 0 1 sq2 -sq2 1 -sq2 sq3;
                           -sq2 1 0 1 sq3 sq2 1 -sq2;
                           1 sq2 1 0 -sq2 -sq3 sq2 1;
                           1 -sq2 sq3 -sq2 0 1 sq2 1;
                           sq2 1 sq2 -sq3 1 0 1 -sq2;
                           -sq3 -sq2 1 sq2 sq2 1 0 1;
                           sq2 sq3 -sq2 1 1 -sq2 1 0];
                QBp = 10^1/nAg; % original 10^1/nAg;
                QBdp = 10^0/nAg; % original 10^0/nAg;
                ptilde_p = QBp; % terminal QBp  original 10/nAg;
                ptilde_v = QBdp; % terminal QBdp original 1/nAg;
                QB = [QBp*eye(DIM) zeros(DIM);
                      zeros(DIM)   QBdp*eye(DIM)];
                R = 1*eye(nAg*DIM);
                kF = 2;
                kA = 0.25;
                kr = 250;
                ka = 100;
                T = 20;
                T0 = (0:dt:T)';
                x0p = 3*[-2 1 8  -3 -1 1  2 -2 5  5 1 5  ...
                        -1 -6 2  2 2 2  -1 -5 -1  5 2 3];   
                x0v = 5*[2*vel vel 0  2*vel 0 0  0 0 0  3*vel -vel vel  ...
                         0 0 vel  vel vel vel  0 0 0 -2*vel -5*vel 2*vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel vel vel];
                
                % ------------------------------------------------
%                 x0v = 5*[2*vel vel 0  2*vel 0 0  0 0 0  3*vel -vel vel  ...
%                          0 0 vel  vel vel vel  0 0 0 -2*vel -5*vel 2*vel];
%                 % helix curve 
%                 rr = 15;
%                 cc = 2*vel;
%                 x_hel = rr*cos(T0);
%                 y_hel = rr*sin(T0);
%                 z_hel = cc*T0;
%                 dx_hel = -y_hel;
%                 dy_hel = x_hel;
%                 dz_hel = cc*ones(length(T0),1);
%                 for t_p = 1:length(T0)
%                     p_hel = Rzx*[x_hel(t_p) y_hel(t_p) z_hel(t_p)]';
%                     x_hel(t_p) = p_hel(1);
%                     y_hel(t_p) = p_hel(2);
%                     z_hel(t_p) = p_hel(3);
%                     v_hel = Rzx*[dx_hel(t_p) dy_hel(t_p) dz_hel(t_p)]';
%                     dx_hel(t_p) = v_hel(1);
%                     dy_hel(t_p) = v_hel(2);
%                     dz_hel(t_p) = v_hel(3);
%                 end
%                 Xdesp = repmat(xd0p,length(T0),1) + [x_hel y_hel z_hel];
%                 Xdesv = [dx_hel dy_hel dz_hel];
                x0v = 5*[2*vel vel 0  2*vel 0 0  0 0 0  3*vel -vel vel  ...
                         0 0 vel  vel vel vel  0 0 0 -2*vel -5*vel 2*vel];
                % tanh curve 
                rr = 10;
                cc = 2*vel;
                kc = 10;
                x_hel = cc*T0;
                y_hel = rr*tanh(kc*(T0-T/2));
                z_hel = 0*T0;
                dx_hel = cc*ones(length(T0),1);
                dy_hel = (kc*rr/(cosh(T0-T/2)).^2)';
                dz_hel = 0*T0;
                for t_p = 1:length(T0)
                    p_hel = Rzy*[x_hel(t_p) y_hel(t_p) z_hel(t_p)]';
                    x_hel(t_p) = p_hel(1);
                    y_hel(t_p) = p_hel(2);
                    z_hel(t_p) = p_hel(3);
                    v_hel = Rzy*[dx_hel(t_p) dy_hel(t_p) dz_hel(t_p)]';
                    dx_hel(t_p) = v_hel(1);
                    dy_hel(t_p) = v_hel(2);
                    dz_hel(t_p) = v_hel(3);
                end
                Xdesp = repmat(xd0p,length(T0),1) + [x_hel y_hel z_hel];
                Xdesv = [dx_hel dy_hel dz_hel];
                
                % I've commented the 2 lines below to use the helix/ tanh
                % ------------------------------------------------
                %Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                %Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                T0e = T0(end);
                x11p = [];
                phi = pi/4;
                th = pi/2;
%                 for i = 1:nAg/2
%                     x11p = [x11p...
%                         vel*T0e + sq3/2*dd*cos(phi)*cos((i-1)*th) ...
%                         vel*T0e + sq3/2*dd*cos(phi)*sin((i-1)*th)...
%                         vel*T0e + sq3/2*dd*sin(phi)...
%                         vel*T0e + sq3/2*dd*cos(-phi)*cos((i-1)*th) ...
%                         vel*T0e + sq3/2*dd*cos(-phi)*sin((i-1)*th)...
%                         vel*T0e + sq3/2*dd*sin(-phi)];
%                 end
                %x11v = vel*ones(1,DIM*nAg); % straight line
                pcfin = Rzy*[cc*T0e rr*tanh(kc*T0e/2) 0]'; % tanh curve
                % pcfin = Rzx*[rr*cos(T0e) rr*sin(T0e) cc*T0e]';
                for i = 1:nAg/2
                    x11p = [x11p...
                        pcfin(1) + sq3/2*dd*cos(phi)*cos((i-1)*th) ...
                        pcfin(2) + sq3/2*dd*cos(phi)*sin((i-1)*th)...
                        pcfin(3) + sq3/2*dd*sin(phi)...
                        pcfin(1) + sq3/2*dd*cos(-phi)*cos((i-1)*th) ...
                        pcfin(2) + sq3/2*dd*cos(-phi)*sin((i-1)*th)...
                        pcfin(3) + sq3/2*dd*sin(-phi)];
                end
                x11v = kron(ones(nAg,1),[cc/2 cc/2 cc/sqrt(2)]')';       % tanh
                % x11v = kron(ones(nAg,1),Rzx*[-rr*sin(T0e) rr*cos(T0e) cc]')';
                x11 = [ x11p x11v ];
                kT = 5;
                
            case 2
                % full cube
                sq2 = sqrt(2);
                sq3 = sqrt(3);
                dijs = dd*[0 1 sq2 1 1 sq2 sq3 sq2;
                           1 0 1 sq2 sq2 1 sq2 sq3;
                           sq2 1 0 1 sq3 sq2 1 sq2;
                           1 sq2 1 0 sq2 sq3 sq2 1;
                           1 sq2 sq3 sq2 0 1 sq2 1;
                           sq2 1 sq2 sq3 1 0 1 sq2;
                           sq3 sq2 1 sq2 sq2 1 0 1;
                           sq2 sq3 sq2 1 1 sq2 1 0];
                QBp = 10^1/nAg;
                QBdp = 10^0/nAg;
                ptilde_p = 10/nAg; % terminal QBp 
                ptilde_v = 1/nAg; % terminal QBdp
                QB = [QBp*eye(DIM) zeros(DIM);
                      zeros(DIM)   QBdp*eye(DIM)];
                R = 1*eye(nAg*DIM);
                kF = 2;
                kA = 0.25;
                kr = 250;
                ka = 100;
                T = 20;
                T0 = (0:dt:T)';
                x0p = 3*[-2 1 8  -3 -1 1  2 -2 5  5 1 5  ...
                        -1 -6 2  2 2 2  -1 -5 -1  5 2 3];   
                x0v = 5*[2*vel vel 0  2*vel 0 0  0 0 0  3*vel -vel vel  ...
                         0 0 vel  vel vel vel  0 0 0 -2*vel -5*vel 2*vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel vel vel];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                T0e = T0(end);
                x11p = [];
                phi = pi/4;
                th = pi/2;
                for i = 1:nAg/2
                    x11p = [x11p...
                        vel*T0e + sq3/2*dd*cos(phi)*cos((i-1)*th) ...
                        vel*T0e + sq3/2*dd*cos(phi)*sin((i-1)*th)...
                        vel*T0e + sq3/2*dd*sin(phi)...
                        vel*T0e + sq3/2*dd*cos(-phi)*cos((i-1)*th) ...
                        vel*T0e + sq3/2*dd*cos(-phi)*sin((i-1)*th)...
                        vel*T0e + sq3/2*dd*sin(-phi)];
                end
                x11v = vel*ones(1,DIM*nAg);
                x11 = [ x11p x11v ];
                kT = 5;
        end
    
end

for i = 1:nAg
    for j = 1:i-1
        if dijs(i,j) > 0
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end

%% desired trajectories

% just in case ...
TT = length(T0);
Tdes = T0;

% desired trajectory

Udes = repmat(0,length(T0),NI);

% initial trajectory

U0 = repmat(0,length(T0),NI);

%% *** CENTRALIZED CONTROLLER DEVISED WITH PRONTO ***

%% PD controller
wn = 3;
zeta = 0.7;
kv = 2*zeta*wn;
kp = wn^2;

% Kpdt = [ kp*eye(NI) , kv*eye(NI) ]';
% Kpd_ = Kpdt(:)';

Kpd = [ kp*eye(NI) , kv*eye(NI) ];
Kpd_ = reshape(Kpd',1,2*NI*NI);

Kr = repmat(Kpd_,length(T0),1);

%% PRONTO setup

format short g


Wt = [];
Wlt  = Xdes; 


% % % NS2 = NS*NS;
NN12 = NS*(NS+1)/2;

% options for 'sim'
opts = simset('Solver','ode45', ...
	      'AbsTol',1e-8,'RelTol',1e-6, ...
              'Trace','minstep');  % to catch errors ... later
          
% turn off some warnings that are not applicable to our MDLs
if 0
  % load_system('LQR_KrS');
  load_system('nonl_KS');
  load_system('Prq1_KvS');
  load_system('Prq2_KvS');
  load_system('linK_zvS');
  % set_param('LQR_KrS',  'InitInArrayFormatMsg', 'None');
  set_param('nonl_KS',  'InitInArrayFormatMsg', 'None');
  set_param('Prq1_KvS', 'InitInArrayFormatMsg', 'None');
  set_param('Prq2_KvS', 'InitInArrayFormatMsg', 'None');
  set_param('linK_zvS', 'InitInArrayFormatMsg', 'None');
  if 0
    % save_system('LQR_KrS');
    save_system('nonl_KS');
    save_system('Prq1_KvS');
    save_system('Prq2_KvS');
    save_system('linK_zvS');
  end
end
         
fig_base = 1000;


% get screen size for plotting
scrsz = get(0,'ScreenSize');

wscr = scrsz(3);
hscr = scrsz(4);

% figure out dimensions of a default figure
%
%   by experimentation on a MacBook Pro
%   (manual description is bad, with ambiguous info to be found!)
%     Position describes the 'drawable area'
%     OuterPosition describes the overall window of the figure
%   Mac MATLAB doesn't include any 'outer frame' as far as I can tell
%
h = figure; drawnow
pos  = get(h,'Position');       % of automatically opened figure
opos = get(h,'OuterPosition');  % includes space for frame & menu and tool bars
set(h,'toolbar','none','menubar','none'); drawnow
ontpos = get(h,'OuterPosition'); % how much smaller?
ntpos = pos;  ntpos(2) = ntpos(2) + opos(4) - ontpos(4);
% set(h,'Position',ntpos,'toolbar','none','menubar','none'); drawnow
close(h);
% the drawnow commands seem to be necessary to set the position parameters!!!

wfig = pos(3);
hfig = pos(4);
hdfig = opos(4) -  pos(4);    % size of header including menu and tool bars
tpfig = opos(2) + opos(4);    % location of top of figure
hdntfig = ontpos(4) - pos(4); % header size without tool bars

% make X,U and Z,V figs have same aspect ratio
wnfig = round( (wscr - wfig)/2 - 5 );
hnfig = round( hfig/wfig*wnfig );
ytnfig = tpfig - (hnfig + hdntfig);
ybnfig = ytnfig - (hnfig + hdntfig) - 5;

FnPos = ntpos;
XiPos = [0 ytnfig wnfig hnfig];
UiPos = [0 ybnfig wnfig hnfig];
ZiPos = [wscr-wnfig ytnfig wnfig hnfig];
ViPos = [wscr-wnfig ybnfig wnfig hnfig];

% place Descent plot below cost function plot
DescPos = [FnPos(1)+FnPos(3)-wnfig FnPos(2)-hnfig-5 wnfig hnfig];

%% final cost terms initialization
% organize the lower triangle of Pf in a flat manner
C = zeros(2*DIM,NS);
I_DIM = eye(DIM);
for ii = 1:nAg
    C(1:DIM,(ii-1)*DIM+1:(ii-1)*DIM+DIM) = I_DIM;
    C(DIM+1:2*DIM,(nAg+ii-1)*DIM+1:(nAg+ii-1)*DIM+DIM) = I_DIM;
end
C = C/nAg;
Ptilde_p = ptilde_p*eye(DIM);
Ptilde_v = ptilde_v*eye(DIM);
Ptilde = [Ptilde_p zeros(DIM,DIM); zeros(DIM,DIM) Ptilde_v];
Pf0 = nAg*C'*Ptilde*C; % zeros(NS);

%% initial conditions & PRONTO parameters

% grab x00 & x11
x00 = X0(1,:);
% x11 = X0(end,:);

% always start nonl_K at the *same* init cond
ip_opts = simset(opts, 'InitialState', [x00 0]);

Xi = X0;
Ui = U0;

% algorithm params
gam_alf = 0.4;	gam_beta = 0.7;
% MaxIter = 20;	 DescentTol = 1e-4;
% MaxIter = 40;	 DescentTol = 1e-6;
% MaxIter = 200; DescentTol = 1e-6;
% MaxIter = 4;	 DescentTol = 1e-6;
% MaxIter = 30;	DescentTol = 1e-8;
% MaxIter = 40;	DescentTol = 1e-8;
MaxIter = 30;	DescentTol = 1e-8;
% MaxIter = 50;	DescentTol = 1e-8;

max_gam = 1;

Colors = 'bgrcmyk';	nColors = 7;	% for plotting

verbose = 3;	% 0, 1, 2, 3 possible  (also 4---causes a pause
                % each iter)

%% to store results...
Descent = [];
Quad = [];
Vals = [];
MinCost = [];
MaxZV = [];
Method = [];
Gammas = [];
iGammas = [];
CompTimes = [];
gCost_vals = [];
gCost_tr_vals = [];
gCost_in_vals = [];
gCost_fo_vals = [];
XBi_norm = zeros(T,1);   % nAg*/2*|| xB - xB,des ||_QB ^2 for all t in [0,T]
Ui_norm = zeros(T,1);    % 1/2*|| u ||_R ^2 for all t in [0,T]
FOi = zeros(T,1);        % kF/2*(sum_1^n sum_i\neq j sigma(sij))

%% main cycle for PRONTO

for ind = 0:MaxIter
% project and evaluate
  [Ti, aXi, aYi] = ...
    sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Wt Wlt]);

  if 1
    max_dXU = max(abs( [ aXi(:,1:NS)-Xi, aYi(:,1:NI)-Ui ] ));
    fprintf('max_dXU:');
    fprintf(' %g', max_dXU);
    fprintf('\n');
  end

  Xi = aXi(:,1:NS);
  Ui = aYi(:,1:NI);
  if NO > 0
    Yi = aYi(:,NI+1:NI+NO);
  end

  Xf = Xi(end,:);
  dXf = Xf - x11;
  
  Pf = Pf0 + final_FO(Xf,2,DIM,nAg,kF,kA,dijs,kr,ka);
  Pf_ = zeros(1,NS*(NS+1)/2);
  if 1
    ii_ = 0;
      for ii=1:NS
        Pf_(ii_+1:ii_+ii) = Pf(ii,1:ii);
        ii_ = ii_+ii;
      end
  end
  
  
  val = aXi(end,NS+1) + 1/2*dXf*Pf0*dXf' + final_FO(Xf,0,DIM,nAg,kF,kA,dijs,kr,ka);
  %  Ji = val - aXi(:,NS+1);   % cost to go

  if verbose > 1
    Xdes = Wlt(:,1:NWL);
% % %     Udes = Wlt(:,NS+1:NS+NI);

    XB = zeros(length(Xi(:,1)), DIM);
    for kkk = 1:DIM
      for kk = 1:nAg
        XB(:,kkk) = XB(:,kkk) + Xi(:,(kk-1)*DIM+kkk)/nAg;
      end
    end
    
    for kkk = 1:DIM
      figure(fig_base + kkk), set(gcf,'Position',XiPos,'toolbar','none','menubar','none');
      plot(T0,Xdes(:,kkk),'-.',T0,XB(:,kkk),'-')
      grid on, zoom on
      kkk_str = num2str(kkk);
      title(strjoin(strcat('X_B_d_e_s_',kkk_str,{' '},'and',...
          {' '},'X_B_i_',kkk_str)))
    end
    
  end

% keyboard
% pause

  %************************************************
  % get descent direction  zeta = (z,v) 
  %   and compute  Dg(xi).zeta          (= descent)
  %         and    D2g(xi).(zeta,zeta)  (= quad)
  %************************************************

  compTime = cputime;

  rf = Pf0*dXf' + final_FO(Xf,1,DIM,nAg,kF,kA,dijs,kr,ka);
  Prq0 = [Pf_ rf' rf'];
  Popts = simset(opts, 'InitialState', Prq0);

  % heuristic for using FIRST or SECOND order method
  % if ind<3 | descent_toggle < 1,
  if 0 %ind<2      %   | descent_toggle < 1,

    Tb = T0(end) - flipud(T0);
    [Tb,Prq,Kvb] = ...
        sim('Prq1_KvS',Tb,Popts, ...
            [Tb flipud([Xi Ui Kr Wt Wlt])]);
    Method=[Method;1];

  else

    % try using second order method
    try
      Prq_cpu = cputime;
      Tb = T0(end) - flipud(T0);
      [Tb,Prq,Kvb]=sim('Prq2_KvS',Tb,Popts, ...
          [Tb flipud([Xi Ui Kr Wt Wlt])]);
      Prq_cpu = cputime - Prq_cpu;
      fprintf('Prq2_cpu = %g\n',Prq_cpu);

      Method = [Method;2];


    catch last_err
      disp(last_err);
      disp('switching to first order descent');
      Prq_cpu = cputime;
      Tb = T0(end) - flipud(T0);
      [Tb,Prq,Kvb]=sim('Prq1_KvS',Tb,Popts, ...
          [Tb flipud([Xi Ui Kr Wt Wlt])]);
      Prq_cpu = cputime - Prq_cpu;
      fprintf('Prq1_cpu = %g\n',Prq_cpu);

      Method=[Method;1];

    end

  end

  Ki = flipud(Kvb(:,1:(NI*NS)));
  Vo = flipud(Kvb(:,((NI*NS)+1):((NI*NS)+NI)));
  Ri = flipud(Prq(:,(NN12+1):(NN12+NS)));
  Qi = flipud(Prq(:,(NN12+NS+1):(NN12+NS+NS)));
  Pi = flipud(Prq(:,1:NN12));

  [Tb,aZi,aVi]=sim('linK_zvS',T0,opts, ...
                  [T0 Xi Ui Ki Vo Qi Wt Wlt]);
  Zi = aZi(:,1:NS);
  Vi = aVi(:,1:NI);

  Zf = Zi(end,:);
  descent = aZi(end,NS+1) + rf'*Zf';
  quad = aZi(end,NS+2) + Zf*Pf*Zf';

  fprintf('ind: %d: %g\n', ind, val);
  fprintf('  descent: %g   term: %g\n', descent, rf'*Zf');
  fprintf('  quad: %g   term: %g\n', quad, 1/2*Zf*Pf*Zf');


  Vals = [Vals; val];
  Descent = [Descent; descent];
  Quad = [Quad; quad];
  

  % restrict the step length (max_gam)
  %   to some max change in phi, theta (say, 30 degrees) ***special***
  %maxZi = max(max(abs(Zi(:,1))));

  % MinCost = [MinCost; min_cost];
  %MaxZV = [MaxZV; [ maxZi, max(abs(Vi))] ];

  if verbose > 0
    figure(fig_base + kkk+1), set(gcf,'Position',ZiPos,'toolbar','none','menubar','none');
      plot(Ti, Zi)
      grid on, zoom on
     title('Zi')

    figure(fig_base + kkk+2), set(gcf,'Position',ViPos,'toolbar','none','menubar','none');
      plot(Ti, Vi)
      grid on, zoom on
    title('Vi')

    if ind>0
     figure(fig_base + kkk+3), set(gcf,'Position',DescPos,'toolbar','none','menubar','none');
      plot(0:(length(Descent)-1), log10(-Descent),'g-', ...
           0:(length(Descent)-1), log10(-Descent),'r+', ...
           0:(length(Vals)-2), log10(-diff(Vals)),'c-o')
      grid on, zoom on
     title('log10 -Descent')
    end

  end

  if verbose > 3
    pause
  end

  if ind > 0
    % save trace of computation
    eval(sprintf('T%i = T0;',ind));
    eval(sprintf('X%i = Xi;',ind));
    eval(sprintf('U%i = Ui;',ind));
  end

  


  %************************************************
  % backtracking (armijo) line search
  %************************************************
  %
  %   using  gam_alf in ( 0, 1/2 ),  gam_beta in ( 0, 1 )
  % 
  %   accept gamma = gam_beta ^ k, k = 0, 1, ...
  %
  %     g(xi + gamma*zeta) < g(xi) +  gam_alf * descent * gamma
  %

  gam_accept = 0;

  % ***special***
% % %   max_gam = 0.5/maxZi;     % gam=1 corresponds to 0.5 rads in phi, theta max
% % %   max_gam = 1.0/maxZi;     % gam=1 corresponds to 1.0 rads in phi, theta max
  max_gam = 1.0; %***************
  gammai = min(1.0,max_gam);
  % gammai = 1.0;
  gamsi = [];
  gCost = [];
  gCost_tr = [];
  gCost_in = [];
  gCost_fo = [];
  residual_terminal_cost_tr = 0;
  residual_terminal_cost_fo = 0;
  for k = 0:12
    [Ti,aXi,aYi] = ...
      sim('nonl_KS',T0,ip_opts, ...
          [T0 ([Xi Ui]+gammai*[Zi Vi]) Kr Wt Wlt]);
    g_dXf = aXi(end,1:NS) - x11;
    gamsi = [gamsi; gammai];
    % ---------------------------------------------------------------------------
    g_FO_gammai = final_FO(aXi(end,1:NS),0,DIM,nAg,kF,kA,dijs,kr,ka);
    residual_terminal_cost_tr = 1/2*g_dXf*Pf0*g_dXf';
    residual_terminal_cost_fo = g_FO_gammai;
    gCost = [gCost; aXi(end,NS+1) + residual_terminal_cost_tr + residual_terminal_cost_fo];

    disp(sprintf('gCost: %g', gCost(end)));

    if gCost(end) < val  +  gam_alf * descent * gammai
      gam_accept = 1;
      % igam = find(gCost == min(gCost));
      igam = find(gCost == gCost(end));  % silly quick fix
      iGammas = [iGammas; igam];
      gammai = gamsi(igam);
      Gammas = [Gammas; gammai];
      break;
    end

    gammai = gam_beta*gammai;
  end
  
  gCost_vals = [gCost_vals gCost(end)];
  [gCost_tri, gCost_ini, gCost_foi, XBi_norm, Ui_norm, FOi] =...
      computeTrInFoCosts(aXi(:,1:NS),aYi(:,1:NI),Xdes,QB,R,nAg,DIM,T0,...
      dijs,kF,kA,kr,ka);
  gCost_tr_vals = [gCost_tr_vals gCost_tri+residual_terminal_cost_tr];
  gCost_in_vals = [gCost_in_vals gCost_ini];
  gCost_fo_vals = [gCost_fo_vals gCost_foi+residual_terminal_cost_fo];
  
  % check for the costs 
%   crt = gCost_tri+residual_terminal_cost_tr
%   cin = gCost_ini
%   cfo = gCost_foi+residual_terminal_cost_fo
%   crt + cin + cfo
%   gCost(end)
%   pause
  
  if descent > -DescentTol
    disp(sprintf('convergence achieved: descent: %g\n', descent));
    break;
  end

  compTime = cputime - compTime;
  CompTimes = [CompTimes; compTime];

  %*******************************

  if 0 %verbose
    % gams = (0:.1:2.5)*min(1.0,max_gam);
    gams = (0:.1:1.2) * min(1.0,max_gam);
    Cost = zeros(size(gams));
  
    for k = 1:length(gams)
      gamma = gams(k);
      [Ti,aXi,aYi] = ...
        sim('nonl_KS',T0,ip_opts, ...
            [T0 ([Xi Ui]+gamma*[Zi Vi]) Kr Wt Wlt]);
      g_dXf = aXi(end,1:NS) - x11;
      Cost(k) = aXi(end,NS+1) + 1/2*g_dXf*Pf*g_dXf';
    end
  
    figure(fig_base + 10 + ind)
      set(gcf,'Position',FnPos,'toolbar','none','menubar','none');
      plot(gams,min(Cost,2*Cost(1)), ...
  	gamsi, min(gCost, 2*Cost(1)), 'x', ...
  	... % [0 .5], Cost(1)+descent*[0 .5], ...
    	... % [0 1.3], Cost(1)+descent/2*[0 1.3], ...
  	gams, Cost(1)+descent*gams, ...
    	gams, Cost(1)+descent/2*gams, ...
    	gams, Cost(1) + descent*gams + quad*gams.*gams/2, ...
    	gams, Cost(1)+0.4*descent*gams ...
	... % [0 1.3], Cost(1)+0.4*descent*[0 1.3] ...
          )
      title(sprintf('order %d: [X%d U%d] + gamma*[Z%d V%d]',Method(end),ind,ind,ind,ind));
      grid on, zoom on
    drawnow

    %%% jh, 9 may 07
% % %       if 0*1   %% ind<1
% % %         pause
% % %       end
    %%%

    gammaS = -descent/quad;
    disp(sprintf('  gammaS:  %g   predict:  %g', gammaS, Cost(1)-descent^2/quad/2));
  end

  if ~ gam_accept 
    % error( ) ??
    disp('*** unable to satisfy step size --- aborting');
    break;
  end

  disp(sprintf('  gammai = %g', gammai));
  
  % **** history plot of the cost descent ***
  
  if 1
  end


  %************************************************
  % update
  %************************************************

  [Ti,aXi,aYi] = ...
    sim('nonl_KS',T0,ip_opts, ...
        [T0 ([Xi Ui]+gammai*[Zi Vi]) Kr Wt Wlt]);

  Xf = aXi(end,1:NS);
  g_dXf = aXi(end,1:NS) - x11;
  gu_FO_val = final_FO(Xf,0,DIM,nAg,kF,kA,dijs,kr,ka);
  term_val_u = 1/2*g_dXf*Pf0*g_dXf' + gu_FO_val;
  val = aXi(end,NS+1) + term_val_u;

  disp(sprintf('val: %g   final_cost: %g   norm(dXf): %g', val, term_val_u, norm(g_dXf)));
  disp(sprintf('Xf: %g %g', Xf(1:NS)));

  Xi = aXi(:,1:NS);
  Ui = aYi(:,1:NI);
  Ji = val - aXi(:,NS+1);

  Xf = Xi(end,:);
  dXf = Xf - x11;

% % %   if 0    % verbose>1
% % %     figure(fig_base + 1)
% % %       plot(To,Xo(:,[1]),'-.',T0,X0(:,[1]),'-.',Ti,Xi(:,[1]),'-')
% % %       grid on, zoom on
% % %     title('Xo and Xi')
% % %  
% % %     figure(fig_base + 2)
% % %       plot(To,Uo,'-.',T0,U0,'-.',Ti,Ui,'-')
% % %       grid on, zoom on
% % %     title('Uo and Ui')
% % %   end

  % pause % ****

end

disp(sprintf('newt_pdbot CompTime: %g', sum(CompTimes)));


Xopt = Xi;
Uopt = Ui;
Kopt = Ki;

%% graphics
ftsz = 35;                % fontsize
lw = 3;                   % linewidth

grp.XBi_norm = XBi_norm;
grp.Ui_norm = Ui_norm;
grp.XBdes = Xdes;
grp.Xi = Xi;
grp.kT = kT;
grp.gCost_tr_vals = gCost_tr_vals;
grp.gCost_in_vals = gCost_in_vals;
grp.gCost_fo_vals = gCost_fo_vals;
grp.gCost_vals = gCost_vals;
grp.Descent = Descent;
grp.FOi = FOi;
grp.Ui = Ui;
grp.T0 = T0;

% centralized graphics
graphics_main(grp,'centralized',nAg,DIM,ftsz,lw,dijs,R,kr,ka,A)

error('\nFINE PRONTO\n')
%pause


%% *** DISTRIBUTED CONTROL BASED ON VARIATIONAL CALCULUS ***

close all

R = 1*10^0*eye(nAg*DIM); % 
QBp = 1*10^1*eye(DIM)/nAg; % 1*10^1*eye(DIM)/nAg 
QBdp = 1*10^0*eye(DIM)/nAg; % 
QB = [QBp           zeros(DIM);
      zeros(DIM)    QBdp];
kF = 2*10^0; % 
kA = 0.25*10^0; % 
kr = 250*10^0; % 250
ka = 100*10^0; %
alpha = [1 1]; % not used

% gains for the PD controller (pentagon side 5 diagonal 5*2*sin(54�))
% kP_tr = 0.4; 
% kP_fo = 2*10^0; 
% kD_tr = 2*10^0; 
% kD_fo = 2*10^-1; % 2 pent, 1 fullpent
% kD_vel = 1; %8*10^-1;
% gains = [kP_tr kP_fo kD_tr kD_fo kD_vel]/1; 

% gains for the PD controller (cube)
kP_tr = 3*10^-1*nAg; % 1*10^-1*nAg
kD_tr = 3*10^0*nAg; % 1.4*10^0*nAg    = (kD_tr1+kD_tr2/10)
kP_fo = 1.3*10^0; % 1.3*10^0
kD_fo = 1*10^0; % 1*10^0
kP_al = 3*10^-1; % 3*10^-1
gains = [kP_tr kD_tr kP_fo kD_fo kP_al]/1; 

gammaOBS_P = 180;%200*max(eig(QBp))*kP_tr; % 200
gammaOBS_D = 170;%150*max(eig(QBdp))*kD_tr; % 210


% framework definition
dist_coeff = 2;

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
% F = eye(nAg)*eta_star+(D^-1*A)*(1-eta_star);
F = get_doubly_stoch(A);

% integration of the differential equation representing the dynamics

% parameters
params.kr = kr;
params.ka = ka;
params.dt = dt;
params.nAg = nAg;
params.DIM = DIM;
params.QBp = QBp;
params.QBdp = QBdp;
params.R = R;
params.kF = kF;
params.kA = kA;
params.qA = A;
params.gains = gains;
params.Xdes = Xdes;
params.dijs = dijs;
params.G = G;
params.tl = T;
params.F = kron(F,eye(DIM));
params.dist_coeff = dist_coeff;
params.max_dist = max_dijs;
params.T0 = T0;
params.alpha = alpha;
params.gammaOBS_P = gammaOBS_P;
params.gammaOBS_D = gammaOBS_D;

% flags to show completion
flag33 = 0;
flag50 = 0;
flag90 = 0;

% dynamics integration
M = DIM;
x0p_hat = zeros(1,nAg*NI);
x0v_hat = zeros(1,nAg*NI);
for i = 1:nAg
    Ni = neighbors(G,i);
    intval_i = (i-1)*M+1:(i-1)*M+M;
    intval_ii = (i-1)*NI+(i-1)*M+1:(i-1)*NI+(i-1)*M+M;
    x0p_hat(intval_ii) = x0p(intval_i);
    x0v_hat(intval_ii) = x0v(intval_i);
    wii = 1+length(Ni);
    p0hati = wii*x0p_hat(intval_ii);
    v0hati = wii*x0v_hat(intval_ii);
    deni = wii;
    for j_ = 1:length(Ni)
        j = Ni(j_);
        intval_j = (j-1)*M+1:(j-1)*M+M;
        intval_ij = (i-1)*NI+(j-1)*M+1:(i-1)*NI+(j-1)*M+M;
        x0p_hat(intval_ij) = x0p(intval_j);
        x0v_hat(intval_ij) = x0v(intval_j);
        wij = 1+length(neighbors(G,j));
        p0hati = p0hati + wij*x0p_hat(intval_ij);
        v0hati = v0hati + wij*x0v_hat(intval_ij);
        deni = deni + wij;
    end
    p0hati = p0hati/deni;
    v0hati = v0hati/deni;
    Nicompl = setdiff((1:nAg),union(Ni,i));
    for j_ = 1:length(Nicompl)
        j = Nicompl(j_);
        intval_ij = (i-1)*NI+(j-1)*M+1:(i-1)*NI+(j-1)*M+M;
        x0p_hat(intval_ij) = p0hati;
        x0v_hat(intval_ij) = v0hati;
    end
end

[t,xx_history] = ode45(@(t,xx)distr_dyn(t,xx,params),T0',[x0p x0v x0p_hat x0v_hat]');
x_history = xx_history(:,1:NS)
u_history = diff(x_history(:,NI+1:NS))/dt;
u_history = [u_history; u_history(end,:)];

% fase diagram for positions
[gCost_tri, gCost_ini, gCost_foi, XBi_norm, Ui_norm, FOi] =...
      computeTrInFoCosts(x_history(:,1:NS),u_history(:,1:NI),Xdes,...
      QB,R,nAg,DIM,T0,dijs,kF,kA,kr,ka);
  
% parameters
grp.XBi_norm = XBi_norm;
grp.XBdes = Xdes;
grp.Xi = x_history;
grp.kT = kT;
grp.gCost_tr_vals = gCost_tri;
grp.gCost_in_vals = gCost_ini;
grp.gCost_fo_vals = gCost_foi;
grp.gCost_vals = gCost_vals;
grp.T0 = T0;
grp.FOi = FOi;
grp.Ui = u_history;

% distributed graphics
graphics_main(grp,'distributed',nAg,DIM,ftsz,lw,dijs,R,kr,ka,A)
%%
graphics_est(G,NI,M,gammaOBS_P,gammaOBS_D,xx_history,u_history,T0,ftsz,lw)
%%
graphics_centr_vs_distr(Ui,u_history,T0,ftsz,lw,NI)



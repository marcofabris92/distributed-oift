% newt_Mform.m
%
%     for Formation Control and Path Following
%
%     find the trajactories that allow agents to reach a formation
%     while their barycenter follows a desired path
%     and their input energy is minimized
%
%  Marco Fabris and John Hauser
%  June 2018
%
%  ...
%  

global NS NI QB R nAg DIM T0 dijs kF kr ka ftsz lw
DIM = 2;
SIM_ECC19 = 1;
kTgif = 100;

% get system sizes: NS, NI, NO, etc.
[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;
nAg = NI/DIM;
dt = 0.01;
vel = 1;
QB = [10*eye(DIM) zeros(DIM);
      zeros(DIM)  1*eye(DIM)];
R = 1*eye(nAg*DIM);
kF = 0.1;
kr = 100;
ka = 1;
dijs = [0 5 5;
        5 0 5;
        5 5 0];
tl = 0;
x11 = zeros(1,NS);
Rot = @(theta) ([cos(theta) -sin(theta); sin(theta) cos(theta)]);

switch DIM
    %% 2D simulations for ECC19
    case 2
        switch SIM_ECC19
            
            case 1 % Validation nAg == 3
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-2 1  -3 -1  2 -2];   % randomizer 1*(2*(rand(1,DIM*nAg)-0.5));
                x0p = x0p/(2*norm(x0p));
                x0v = 5*[0 -vel  0 -vel  0 -vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel 0];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                T0e = T0(end);
                nd = 10;
                x11p = [T0e nd*sqrt(3)/3 T0e-nd/2 -nd*sqrt(3)/6 T0e+nd/2 -nd*sqrt(3)/6];
                paff = [(x11p(1)+x11p(3)+x11p(5))/3 (x11p(2)+x11p(4)+x11p(6))/3]';
                for i = 1:nAg
                    iint = (i-1)*DIM+1:(i-1)*DIM+DIM;
                    x11p(iint) = ((Rot(0*pi/180)*(x11p(iint)'-paff))+paff)';
                end
                x11v = [0 0 0 0 0 0];
                x11 = [ x11p x11v ];
                kT = 5;
                
            case 2 % Complicated desired trajectory nAg == 3
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-0.5 -0.5  0 0  6 6];
                x0v = zeros(1,DIM*nAg);
                xd0p = zeros(1,DIM);
                tanh_T0 = 2*tanh(T0-t1/2); % tanh-like curve
                Xdesp = repmat(xd0p,length(T0),1) + [T0*vel  tanh_T0*vel];
                Xdesv = [vel*ones(length(T0),1)  8*exp(2*T0+2*t1/2)./(exp(2*T0)+exp(t1)).^2];
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 5;
                
            case 3 % Initial conditions in a subspace nAg == 3
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-1 0  0 0  1 0]; 
                x0v = [0 0  0 0  0 0];
                xd0p = zeros(1,DIM);
                xd0v = [vel 0];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 5;
                
            case 4 % Equilibria for polygons nAg = [5 6 8 12 20]
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = zeros(1,DIM*nAg);
                for k = 0:nAg-1
                    x0p(k*DIM+1) = cos(2*k*pi/nAg)+rand;
                    x0p(k*DIM+2) = sin(2*k*pi/nAg)+rand;
                end
                x0v = zeros(1,DIM*nAg);
                xd0p = zeros(1,DIM);
                xd0v = zeros(1,DIM);
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 2;
                
        end
  %% 3D simulations for ECC19
    case 3
        switch SIM_ECC19
            
            case 0 % Validation nAg == 1
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-2 1 0];
                x0p = 12*x0p/norm(x0p);
                x0v = 5*[0 -vel vel];
                xd0p = zeros(1,DIM);
                xd0v = [vel vel vel];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 5;
            
            case 1 % Validation nAg == 4
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-2 1 0  -3 -1 1  2 -2 2  1 3 3];
                x0p = 12*x0p/norm(x0p);
                x0v = 5*[0 -vel vel  0 -vel 0  0 -vel 2*vel  0 0 0];
                xd0p = zeros(1,DIM);
                xd0v = [vel vel vel];
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 1; %5;
                
            case 2 % Complicated desired trajectory nAg == 4
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-5 -5 0  0 0 0  6 6 0  -2 2 0];
                x0p = [-5 -5 0  0 0 2  6 6 0  -2 2 0]; % variant
                x0v = zeros(1,DIM*nAg);
                xd0p = zeros(1,DIM);
                
                % helix curve
                rr = 15;
                cc = 2*vel;
                x_hel = rr*cos(T0);
                y_hel = rr*sin(T0);
                z_hel = cc*T0;
                dx_hel = -y_hel;
                dy_hel = x_hel;
                dz_hel = cc*ones(length(T0),1);
                
                Xdesp = repmat(xd0p,length(T0),1) + [x_hel y_hel z_hel];
                Xdesv = [dx_hel dy_hel dz_hel];
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 5;
                
            case 3 % Initial conditions in a subspace nAg == 4
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = [-1 0 0  0 0 0  1 0 0  2 0 0]; 
                x0v = [0 0 0  0 0 0  0 0 0  0 0 0];
                xd0p = zeros(1,DIM);
                xd0v = [vel vel/200 0];
                Xdesp = repmat(xd0p,length(T0),1) + [T0*xd0v(1) T0.^2*xd0v(2) T0*xd0v(3)];
                Xdesv = [xd0v(1)*ones(length(T0),1)  2*T0*xd0v(2) xd0v(3)*ones(length(T0),1)];
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 5;
            
            case 4 % Equilibria for polygons nAg = [5 6 8 12 20]
                t1 = 20;
                T0 = (0:dt:t1)';
                x0p = zeros(1,DIM*nAg);
                for k = 0:nAg-1
                    x0p(k*DIM+1) = cos(2*k*pi/nAg)+rand;
                    x0p(k*DIM+2) = sin(2*k*pi/nAg)+rand;
                    x0p(k*DIM+3) = rand;
                end
                x0v = zeros(1,DIM*nAg);
                xd0p = zeros(1,DIM);
                xd0v = zeros(1,DIM);
                Xdesp = repmat(xd0p,length(T0),1) + T0*xd0v;
                Xdesv = repmat(xd0v,length(T0),1);
                Xdes = [ Xdesp Xdesv ];
                X0p = repmat(x0p,length(T0),1) + T0*x0v;
                X0v = repmat(x0v,length(T0),1); 
                X0 = [ X0p X0v ];
                kT = 2;    
        end
end

%% straight line trajectory

% x0p = zeros(1,DIM*nAg);
% x0p(1) = -2;
% x0p(5) = 1;
% x0p(7) = -3;

% just in case ...
Tdes = T0;

% desired trajectory

Udes = repmat(0,length(T0),NI);

% initial trajectory

U0 = repmat(0,length(T0),NI);


% load X0U0  % used saved X0 U0


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


%%

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

if 0
  x00 = [repmat(Xdes(1,1:DIM),1,nAg)...
         repmat(Xdes(1,DIM+1:2*DIM),1,nAg)]; % Xdes(1,:);  % cost specific! 
  u00 = Udes(1,:);  % cost specific! 

  % for terminal cost
  x11 = [repmat(Xdes(end,1:DIM),1,nAg)...
         repmat(Xdes(end,DIM+1:2*DIM),1,nAg)]; % Xdes(1,:);  % cost specific!
  u11 = Udes(end,:);  % cost specific!  

  % use, as INITIAL TRAJECTORY, an equilbrium traj
  %X0 = repmat(x00,length(T0),1);
  %U0 = repmat(u00,length(T0),1);

  %%

  % get linearization about initial and final conditions
  [~,~,A0,B0] = dynamics_m(x00,u00);
  [~,~,A1,B1] = dynamics_m(x11,u11);

  % get params for regulator and cost

  % for terminal cost
  [~, ~, ~, Q, ~, R] = cost_m(x11, u11, [Xdes(end,1:DIM) Xdes(end,DIM+1:2*DIM)], true);
  % [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt);

  % get terminal cost P1  (call it Pf and Pf_)
  [Kf,Pf,Ef] = lqr(A1,B1,Q,R);
end

% organize the lower triangle of Pf in a flat manner
C = zeros(2*DIM,NS);
I_DIM = eye(DIM);
for ii = 1:nAg
    C(1:DIM,(ii-1)*DIM+1:(ii-1)*DIM+DIM) = I_DIM;
    C(DIM+1:2*DIM,(nAg+ii-1)*DIM+1:(nAg+ii-1)*DIM+DIM) = I_DIM;
end
C = C/nAg;
ptilde_p = 10;
ptilde_v = 1;
Ptilde_p = ptilde_p*eye(DIM);
Ptilde_v = ptilde_v*eye(DIM);
Ptilde = [Ptilde_p zeros(DIM,DIM); zeros(DIM,DIM) Ptilde_v];
Pf0 = C'*Ptilde*C; % zeros(NS);


% for LQR
% % % if 0
% % % [Qreg, Rreg] = QR_params_m;
% % % end

% % % if 0
% % % % get regulator terminal P
% % % %    (just fix it for desired terminal equil state)
% % % [Kreg,Preg,Ereg] = lqr(A1,B1,Qreg,Rreg);
% % % 
% % % % organize the lower triangle of Preg in a flat manner
% % % Preg_ = [];
% % % for ii=1:NS
% % %   Preg_ = [Preg_ Preg(ii,1:ii)];
% % % end
% % % end
  
  
% the initial trajectory is trivial
% % % Xi = X0;
% % % Ui = U0;

% grab x00 & x11
x00 = X0(1,:);
% x11 = X0(end,:);

% always start nonl_K at the *same* init cond
ip_opts = simset(opts, 'InitialState', [x00 0]);

% regulator also always starts at same init cond --- NOT
% % % Kr_opts = simset(opts, 'InitialState', Preg_);


% in case we are silly and try an open loop initial trajectory
% % % if 0
% % %   Kr = zeros(length(T0),NS*NI);
% % % 
% % %   [Ti, aXi, aYi] = ...
% % %      sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Wt Wlt]);
% % % 
% % %   X0 = aXi(:,1:NS);
% % %   U0 = aYi(:,1:NI);
% % % end

% % % if 0
% % %   figure,plot(T0,X0),grid on, zoom on
% % %   pause
% % % end


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
% MaxIter = 80;	DescentTol = 1e-2;
MaxIter = 50;	DescentTol = 1e-8;

max_gam = 1;

Colors = 'bgrcmyk';	nColors = 7;	% for plotting

verbose = 3;	% 0, 1, 2, 3 possible  (also 4---causes a pause
                % each iter)

% to store results

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
XBi_norm = zeros(tl,1);   % 1/2*|| xB - xB,des ||_QB ^2 for all t in [0,T]
Ui_norm = zeros(tl,1);    % 1/2*|| u ||_R ^2 for all t in [0,T]
FOi = zeros(tl,1);        % kF/2*(sum_1^n sum_i\neq j sigma(sij))

for ind = 0:MaxIter
    
%   [~, ~, ~, Q, ~, R] = cost_m(Xi(1,:), Ui(1,:), [ Xdes(1,:) ]);
%   eig(Q)
%   pause

  %************************************************
  % evaluate the functional (after designing K)
  %************************************************

  % design Kr for Xi,Ui

  % get regulator terminal P
% % %   if 0
% % %   [~,~,Ar,Br] = dynamics_m(Xi(end,:),Ui(end,:));
% % %   [Kreg,Preg,Ereg] = lqr(Ar,Br,Qreg,Rreg);
% % %   Preg_ = [];
% % %   for ii=1:NS
% % %     Preg_ = [Preg_ Preg(ii,1:ii)];
% % %   end
% % %   Kr_opts = simset(opts, 'InitialState', Preg_);
% % %   end

  % integrate Riccati eqn *backward* in time --- note flipud
% % %   if 0
% % %   Tb = T0(end) - flipud(T0);
% % %   [Tb,Pb,Kb] = sim('LQR_KrS',Tb,Kr_opts, [Tb flipud([Xi Ui])]);
% % %   Kr = flipud(Kb(:,1:(NI*NS)));
% % %   end

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
  
  Pf = Pf0 + final_FO(Xf,2);
  Pf_ = zeros(1,NS*(NS+1)/2);
  if 1
    ii_ = 0;
      for ii=1:NS
        Pf_(ii_+1:ii_+ii) = Pf(ii,1:ii);
        ii_ = ii_+ii;
      end
  end
  
  
  val = aXi(end,NS+1) + 1/2*dXf*Pf0*dXf' + final_FO(Xf,0);
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

  rf = Pf0*dXf' + final_FO(Xf,1);
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
  for k = 0:12
    [Ti,aXi,aYi] = ...
      sim('nonl_KS',T0,ip_opts, ...
          [T0 ([Xi Ui]+gammai*[Zi Vi]) Kr Wt Wlt]);
    g_dXf = aXi(end,1:NS) - x11;
    gamsi = [gamsi; gammai];
    % ---------------------------------------------------------------------------
    g_FO_gammai = final_FO(aXi(end,1:NS),0);
    gCost = [gCost; aXi(end,NS+1) + 1/2*g_dXf*Pf0*g_dXf' + g_FO_gammai];

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
      computeTrInFoCosts(aXi(:,1:NS),aYi(:,1:NI),Xdes);
  gCost_tr_vals = [gCost_tr_vals gCost_tri];
  gCost_in_vals = [gCost_in_vals gCost_ini];
  gCost_fo_vals = [gCost_fo_vals gCost_foi];
  
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
  gu_FO_val = final_FO(Xf,0);
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

%%

ftsz = 22;                % fontsize
lw = 3;                   % linewidth

fig_base = fig_base + 1000;

%Xdes = Wlt(:,1:NS);
%Udes = Wlt(:,NS+1:NS+NI);

XBdes = Wlt;
XB = zeros(length(Xi(:,1)), DIM);
for kkk = 1:DIM
  for kk = 1:nAg
    XB(:,kkk) = XB(:,kkk) + Xi(:,(kk-1)*DIM+kkk)/nAg;
  end
end

kTgif = 2;
for ii = 1:kTgif-1
    traj_fig = figure(fig_base + kkk+1); set(gcf, 'Position', XiPos+(kkk+1)*[10 -10 30 5]);

    % displying trajectory tracking
    Trajdes = XBdes;
    Ltraj = length(Xi(:,1));
    color_traj = [0 153 76]/255;
    color_bary = [128 255 0]/255;

    if DIM == 2
      plot(Trajdes(:,1),Trajdes(:,2),'k');
      hold on


      for k = 0:nAg-1
        plot(Xi(:,DIM*k+1),Xi(:,DIM*k+2), 'color', color_traj, 'linewidth', 2)
        for kt = 1:kT
          if kt == 1
            tt = 1; %1+floor(ii/kTgif*Ltraj); 
          else
            tt = round((kt-1)/(kT-1)*Ltraj);
          end
          plot(Xi(tt,DIM*k+1),Xi(tt,DIM*k+2), '*', 'color',...
              [1 0 0]/(kT-kt+1)+([0 0 255]/255)*(1-1/(kT-kt+1)),...
              'linewidth', 3)
          txt = strcat('$',num2str(k+1),'$');
          text(Xi(tt,DIM*k+1),Xi(tt,DIM*k+2),txt,'HorizontalAlignment',...
              'right','interpreter','latex','fontsize',ftsz)
        end
      end

      plot(XB(:,1),XB(:,2), 'color', color_bary, 'linewidth', 1);
    end

    if DIM == 3
      plot3(Trajdes(:,1),Trajdes(:,2),Trajdes(:,3),'k');
      hold on

      for k = 0:nAg-1
        plot3(Xi(:,DIM*k+1),Xi(:,DIM*k+2),Xi(:,DIM*k+3), 'color', color_traj, 'linewidth', 2)
        for kt = 1:kT
          if kt == 1
            tt = 1; %1+floor(ii/kTgif*Ltraj); 
          else
            tt = round((kt-1)/(kT-1)*Ltraj);
          end
          plot3(Xi(tt,DIM*k+1),Xi(tt,DIM*k+2),Xi(tt,DIM*k+3), '*', 'color',...
              [1 0 0]/(kT-kt+1)+([0 0 255]/255)*(1-1/(kT-kt+1)),...
              'linewidth', 3)
          txt = strcat('$',num2str(k+1),'$');
          text(Xi(tt,DIM*k+1),Xi(tt,DIM*k+2),Xi(tt,DIM*k+3),txt,...
              'HorizontalAlignment', 'right','interpreter','latex',...
              'fontsize',ftsz)
        end
      end

      plot3(XB(:,1),XB(:,2),XB(:,3), 'color', color_bary, 'linewidth', 1);
    end

    % showing polygons or polyhedra
    color_ini = [1 0 0];
    color_fin = [0 0 1];
    lambda = linspace(0,1);
    for kt = 1:kT
      if kt == 1
        tt = 1+floor(ii/kTgif*Ltraj); %tt = 1; % <------------------------------
     else
        tt = round((kt-1)/(kT-1)*Ltraj);
     end
     if DIM == 2
       plot(XB(tt,1),XB(tt,2),'+','color',...
           color_ini/(kT-kt+1)+color_fin*(1-1/(kT-kt+1)), 'linewidth', 1.5);
     elseif DIM == 3
       plot3(XB(tt,1),XB(tt,2),XB(tt,3),'+','color',...
           color_ini/(kT-kt+1)+color_fin*(1-1/(kT-kt+1)), 'linewidth', 1.5);
     end
      for i = 0:nAg-1
        for j = 0:i-1
          segment = zeros(DIM,length(lambda));
          for k = 1:length(lambda)
            segment(:,k) = (Xi(tt,i*DIM+1:i*DIM+DIM)*lambda(k)+...
                Xi(tt,j*DIM+1:j*DIM+DIM)*(1-lambda(k)))';
          end
          if DIM == 2
            plot(segment(1,:), segment(2,:), 'color',...
              color_ini/(kT-kt+1)+color_fin*(1-1/(kT-kt+1)))
          elseif DIM == 3
            plot3(segment(1,:), segment(2,:), segment(3,:), 'color',...
              color_ini/(kT-kt+1)+color_fin*(1-1/(kT-kt+1)))
          end
        end
      end
    end
    
    ax = gca;
    ax.FontSize = 20;

    hold off
    grid on, zoom on
    %title('Trajectory tracking + formation achieved')
    axis equal
    set(gca,'fontsize', ftsz)
    xlabel('$x$ [m]','Interpreter',...
        'latex', 'FontSize', ftsz, 'linewidth', lw);
    ylabel('$y$ [m]','Interpreter',...
        'latex', 'FontSize', ftsz, 'linewidth', lw);
    if DIM == 3
        zlabel('$z$ [m]','Interpreter',...
            'latex', 'FontSize', ftsz, 'linewidth', lw);
    end
    set(gca,'TickLabelInterpreter','latex')

    figure(traj_fig);
    axis equal
    filename = strcat('traj-',num2str(ii));
    path_office = 'C:\Users\Marco\Desktop\traj_2D_formation';
    formattype = 'pdf';
    %saveas(traj_fig,fullfile(path_office,filename),formattype)
    
    grid on, zoom on
    axis equal
    set(gca,'fontsize', ftsz)
    xlabel('$x$ [m]','Interpreter',...
        'latex', 'FontSize', ftsz, 'linewidth', lw);
    ylabel('$y$ [m]','Interpreter',...
        'latex', 'FontSize', ftsz, 'linewidth', lw);
    if DIM == 3
        zlabel('$z$ [m]','Interpreter',...
            'latex', 'FontSize', ftsz, 'linewidth', lw);
    end
    set(gca,'TickLabelInterpreter','latex')

end

% printing info about distances
p = Xi(end,1:DIM*nAg)'
pp = zeros(DIM,nAg);
for i = 1:nAg
  pp(:,i) = p((i-1)*DIM+1:i*DIM);
end
dists = zeros((nAg*(nAg-1))/2,1);
ij = 1; 
for i = 1:nAg
  for j = 1:i-1
    dists(ij) = norm(pp(:,i)-pp(:,j));
    ij = ij+1;
  end
end
dists


ax = gca;
ax.FontSize = 20;

hold off
grid on, zoom on
%title('Trajectory tracking + formation achieved')
axis equal
set(gca,'fontsize', ftsz)
xlabel('$x$ [m]','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
ylabel('$y$ [m]','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
if DIM == 3
    zlabel('$z$ [m]','Interpreter',...
        'latex', 'FontSize', ftsz, 'linewidth', lw);
end
set(gca,'TickLabelInterpreter','latex')


%% cost descent
figure(fig_base + kkk+2), set(gcf, 'Position', XiPos+(kkk+2)*[0 -80 30 10]);
pk = plot(0:length(gCost_vals)-1, log10(gCost_vals),'k-o', 'MarkerSize',5,'linewidth',lw/2);
hold on
pb = plot(0:length(gCost_tr_vals)-1, log10(gCost_tr_vals),'b-o', 'MarkerSize',5,'linewidth',lw/2);
pg = plot(0:length(gCost_in_vals)-1, log10(gCost_in_vals),'g-o', 'MarkerSize',5,'linewidth',lw/2);
pr = plot(0:length(gCost_fo_vals)-1, log10(gCost_fo_vals),'r-o', 'MarkerSize',5,'linewidth',lw/2);
grid on
zoom on
set(gca,'fontsize', ftsz)
xlabel('$\#$ iterations','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
ylabel('$\log_{10}(h)$','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
set(gca,'TickLabelInterpreter','latex')
legend([pk pb pg pr],{'$h=h^{tr}+h^{in}+h^{fo}$','$h^{tr}$','$h^{in}$','$h^{fo}$'},...
    'location','northeast','interpreter','latex')

%% cost descent (1st order)

figure(fig_base + kkk+3), set(gcf, 'Position', XiPos+(kkk+2)*[0 -80 30 10]);
plot(0:length(Descent)-1, log10(-Descent),'k-o', 'MarkerSize',5)
grid on
zoom on
set(gca,'fontsize', ftsz)
xlabel('$\#$ iterations','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
ylabel('$\log_{10}(-Dh \cdot \zeta)$','Interpreter',...
    'latex', 'FontSize', ftsz, 'linewidth', lw);
set(gca,'TickLabelInterpreter','latex')
hold off

%% tracking vs formation regimes
figure(fig_base + kkk+4), set(gcf, 'Position', XiPos+(kkk+4)*[0 -80 30 10]);
grid on
hold on
[t_settl,t_settl_2,t_settl_3,h_TR,h_FO] = tr_vs_fo(T0,XBi_norm,...
    gCost_tr_vals(end),FOi,gCost_fo_vals(end),1);
average_tracking_time = -h_TR*T0(end)/(h_FO-h_TR);
average_formation_time = h_FO*T0(end)/(h_FO-h_TR);
fprintf('Average tracking time vs average formation time:\n')
display([average_tracking_time average_formation_time])

%% average & total control for each agent
figure(fig_base + kkk+5), set(gcf, 'Position', XiPos+(kkk+5)*[0 -80 30 10]);
grid on
hold on
average_u(T0,Ui,gCost_in_vals(end),t_settl,t_settl_2,t_settl_3);

%% motion decomposition simulations



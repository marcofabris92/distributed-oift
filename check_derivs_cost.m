% check derivatives of cost function
%
%   MF system
%

% check_derivs_cost definitely needs (x,u) and (z,v) ... and these 
% need to be chosen in a manner that might help understand whether 
% things are working ok
% naturally, you can (& should) choose these wisely ...


% [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt)

% clear all, close all, clc

% x = [  p1x p1y p2x p2y v1x v1y v2x v2y ]
xxx = [   0  0.5  0 -0.5  0   0   0   0 ]';
NS = length(xxx);

% u = [  u1x u1y u2x u2y ];
uuu = [   0   0   0   0 ]';

% wlt   = [0 0 0 0 0 0 ];
wlt = [ 0 0 0 0 0 0 0 0 ];

% x = [  p1x p1y p2x p2y v1x v1y v2x v2y ]
zzz = [  1.0  0  0  0       0   0   0   0 ]';
zz0 = [  1.0  0  0  0       0   0   0   0 ]';
% zzz = rand(4,1)/3

% u = [  u1x u1y u2x u2y ];
vvv = [   0   0   0   0 ]';

eps_s = (-1:.001:1)';
% eps_s = (-.5:1e-4:-.4)';

lxu_s = [];
Dlxu_s = [];
Dlxu_zv_s = [];
D2lxu_zv_s = [];

eigQ_s = [];

for ii = 1:length(eps_s)

  epsi = eps_s(ii);
%   disp(sprintf('eps : %g\n',epsi));


% [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt)
  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(xxx+epsi*zzz,uuu+epsi*vvv,wlt);

  lxu_s = [ lxu_s; lxu ];

  Dlxu = [ lxu_x lxu_u ];
  Dlxu_s = [Dlxu_s; Dlxu];
  Dlxu_zv_s = [ Dlxu_zv_s;  Dlxu*[zzz; vvv] ];

  % Build hessian 
  D2lxu = [ lxu_x_x lxu_x_u; 
            lxu_x_u' lxu_u_u];

  D2lxu_zv_s = [D2lxu_zv_s; (D2lxu*[zzz;vvv])' ];

  eigQ_s = [ eigQ_s; eig(lxu_x_x(1:4,1:4))' ];
  
end

% compute forward difference approx

Dlxu_zv_s_fd = diff(lxu_s)  ./  repmat(diff(eps_s),1,length(lxu)) ;
Dlxu_zv_s_fd = [ Dlxu_zv_s_fd; Dlxu_zv_s_fd(end,:) ];

D2lxu_zv_s_fd = diff(Dlxu_s) ./ repmat(diff(eps_s),1,length(Dlxu)) ;
D2lxu_zv_s_fd = [D2lxu_zv_s_fd; D2lxu_zv_s_fd(end,:) ];


% value
figure
  plot(eps_s, lxu_s)
  grid on, zoom on
title('l(x+\\epsilon z,u+\\epsilon v)')

% first derivative
figure
  hold on
  plot(eps_s, Dlxu_zv_s, 'linewidth', 3 )
  plot(eps_s, Dlxu_zv_s_fd, '--', 'linewidth', 2)
  grid on, zoom on
title(sprintf('d/d\\epsilon l(x+\\epsilon z,u+\\epsilon v)'))


return



for ii =  1:length(xxx)

  figure
    hold on
    plot(eps_s, D2lxu_zv_s(:,ii), 'linewidth', 3)
    plot(eps_s, D2lxu_zv_s_fd(:,ii), '--', 'linewidth', 2)
    grid on, zoom on
  title(sprintf('d/d\\epsilon dl_{x%d}(x+\\epsilon z,u+\\epsilon v)',ii))

end


return


if 0
for ii =  1:length(uuu)

  figure
    plot(eps_s, [ D2lxu_zv_s(:,ii+NS) D2lxu_zv_s_fd(:,ii+NS)])
    grid on, zoom on
  title(sprintf('d/d\\epsilon dl_{u%d}(x+\\epsilon z,u+\\epsilon v)',ii))

end
end


return

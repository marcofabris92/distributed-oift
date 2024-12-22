% compare the two functions sig0 and sig1
%
% JH jul 18 boulder

dd = 2;
% dd = 1;

% first look at the functions themselves, then the derivs
%
% we learn:
%
%   * each function has infinite slope at one of the endpoints
%
%   * scaling the functions looks like a good thing

ss = (0:1e-4:2*dd^2)';

sig0 = sig(ss,dd,0);

sig1 = sig(ss,dd,1);

% scale = 1/dd^2/(sqrt(2)-1)^2;
figure, plot(ss, [ sig0 sig1 ]), grid on, zoom on, title('\sigma_0(s) & \sigma_1(s)')


% now, check the derivatives

% avoid the infinite slope areas
ss0 = (0:1e-4:1.5*dd^2)';
ss1 = (0.5*dd^2:1e-4:2*dd^2)';

if 1  % slightly larger s range
  ss1 = (0.5*dd^2:1e-4:3*dd^2)';
end

if 1  % look at large(r)  s = r^2  case
  ss1 = (0.5*dd^2:1e-4:8*dd^2)';
end

[sig0, dsig0, ddsig0] = sig(ss0,dd,0);

[sig1, dsig1, ddsig1] = sig(ss1,dd,1);

figure, plot(ss0, sig0, ss1, sig1), grid on, zoom on, title('\sigma_0(s) & \sigma_1(s)')

figure, plot(ss0, dsig0, ss1, dsig1), grid on, zoom on, title('\sigma''_0(s) & \sigma''_1(s)')

if 0    % compare with finite difference approximation
  dsig0_fd = diff(sig0)./diff(ss0);
  dsig0_fd = [ dsig0_fd; dsig0_fd(end) ];
  dsig1_fd = diff(sig1)./diff(ss1);
  dsig1_fd = [ dsig1_fd; dsig1_fd(end) ];
  figure, plot(ss0, [dsig0 dsig0_fd], ss1, [dsig1 dsig1_fd]), grid on, zoom on, title('\sigma''_0(s) & \sigma''_1(s)')
end


figure, plot(ss0, ddsig0, ss1, ddsig1), grid on, zoom on, title('\sigma''''_0(s) & \sigma''''_1(s)')


if 0    % compare with finite difference approximation
  ddsig0_fd = diff(dsig0)./diff(ss0);
  ddsig0_fd = [ ddsig0_fd; ddsig0_fd(end) ];
  ddsig1_fd = diff(dsig1)./diff(ss1);
  ddsig1_fd = [ ddsig1_fd; ddsig1_fd(end) ];
  figure, plot(ss0, [ddsig0 ddsig0_fd], ss1, [ddsig1 ddsig1_fd]), grid on, zoom on, title('\sigma''''_0(s) & \sigma''''_1(s)')
end



% check C function
sigc    = 0*ss;
sig_s   = 0*ss;
sig_s_s = 0*ss;
for ii=1:length(ss)
  [sig_,sig_s_,sig_s_s_] = sigma_m(ss(ii),dd);
  sigc(ii)    = sig_;
  sig_s(ii)   = sig_s_;
  sig_s_s(ii) = sig_s_s_;
end



figure, plot(ss0, sig0, ss1, sig1, ss, sigc), grid on, zoom on, title('\sigma_0(s) & \sigma_1(s)')

figure, plot(ss0, dsig0, ss1, dsig1, ss, sig_s), grid on, zoom on, title('\sigma''_0(s) & \sigma''_1(s)')

figure, plot(ss0, ddsig0, ss1, ddsig1, ss, sig_s_s), grid on, zoom on, title('\sigma''''_0(s) & \sigma''''_1(s)')




return

% just 'cut-and-paste' the following

% now, look at   tau(r) = sigma( r^2 )  and derivatives
%
%   use sigma_m (written in C) to evaluate

dd = 2;

ss = (0:1e-4:2*dd^2)';

rr = sqrt(ss);

sigc    = 0*ss;
sig_s   = 0*ss;
sig_s_s = 0*ss;
for ii=1:length(ss)
  [sig_,sig_s_,sig_s_s_] = sigma_m(ss(ii),dd);
  sigc(ii)    = sig_;
  sig_s(ii)   = sig_s_;
  sig_s_s(ii) = sig_s_s_;
end

% tau(r) = sigma( r^2 )  is a function of  r
tau = sigc;

% chain rule

tau_r = sig_s .* (2*rr);

tau_r_r = sig_s_s .* (2*rr).^2   + sig_s .* (2);

tau_r_r_safe = sig_s_s .* (2*rr).^2   + max(sig_s,0) .* (2);

% finite difference approx

tau_r_fd = diff(tau) ./ diff(rr);
tau_r_fd = [ tau_r_fd; tau_r_fd(end) ];

tau_r_r_fd = diff(tau_r) ./ diff(rr);
tau_r_r_fd = [ tau_r_r_fd; tau_r_r_fd(end) ];


figure, plot(rr,tau), grid on, zoom on, title('\tau(r)')

figure, plot(rr,[tau_r tau_r_fd]), grid on, zoom on, title('\tau_r(r)')

figure, plot(rr,[tau_r_r tau_r_r_fd]), grid on, zoom on, title('\tau_r_r(r)')

figure, plot(rr,[tau_r_r tau_r_r_safe]), grid on, zoom on, title('\tau_r_r(r) & safe')

close all, clear all, clc

d = 10;

c1 = 10;
c2 = 1;

pts = 10^4;

s = linspace(0,d^2,pts);
L_s = length(s);
sig1 = zeros(1,L_s);
dsig1 = zeros(1,L_s);
ddsig1 = zeros(1,L_s);

for k = 1:L_s
    sig1(k) = c1*(1-s(k)/d^2)^3;
    dsig1(k) = c1*(-3/d^2)*(1-s(k)/d^2)^2;
    ddsig1(k) = c1*(6/d^4)*(1-s(k)/d^2);
end

s_ = s;

s = linspace(d^2,4*d^2,pts);
L_s = length(s);
sig2 = zeros(1,L_s);
dsig2 = zeros(1,L_s);
ddsig2 = zeros(1,L_s);

for k = 1:L_s
    sig2(k) = c2*(sqrt(s(k))/d-1)^3;
    dsig2(k) = 1.5*c2/d^3*(sqrt(s(k))-d)^2/sqrt(s(k));
    ddsig2(k) = 0.75*c2/d^3*(s(k)-d^2)/sqrt(s(k))^3;
end

figure
plot(s_, sig1, 'r', 'linewidth', 2)
hold on
plot(s, sig2, 'b', 'linewidth', 2)
title('sig')
grid on, zoom on

figure
plot(s_, dsig1, 'r', 'linewidth', 2)
hold on
plot(s, dsig2, 'b', 'linewidth', 2)
title('dsig')
grid on, zoom on

figure
plot(s_, ddsig1, 'r', 'linewidth', 2)
hold on
plot(s, ddsig2, 'b', 'linewidth', 2)
title('ddsig')
grid on, zoom on





function output2 = graphics_centr_vs_distr(...
    uc,ud,T0,ftsz,lw,NI)

du = uc-ud;

K = length(T0);
normdu = zeros(K,1);
normuc = zeros(K,1);
for k = 1:K
    normdu(k) = (norm(du(k,:),2));
end
avnormdu = mean(normdu);
Y = log(normdu);
X = [T0 ones(K,1)];
theta = (X'*X)^-1*X'*Y;
m = theta(1);
q = theta(2);




figure
grid on
hold on
rhop = zeros(K,1);
mrhop = log(normdu(end)/normdu(1))/T0(end);
mrhopout = -mrhop;
qrhop = log(normdu(1));
qrhopout = exp(qrhop);
[snormdu,sk] = sort(normdu);
rhom = zeros(K,1);
mrhom = log(snormdu(2)/snormdu(1))/(T0(sk(2))-T0(sk(1)));
qrhom = log(snormdu(2))-mrhom*T0(sk(2));
mrhomout = -mrhom;
qrhomout = exp(qrhom);
output1 = [mrhopout qrhopout mrhomout qrhomout]';
for k = 1:K
    rhop(k) = exp(mrhop*T0(k)+qrhop);
    rhom(k) = exp(mrhom*T0(k)+qrhom);
end

% h3 = plot(T0,log(rhop),'r','linewidth',lw);
% h4 = plot(T0,log(rhom),'b','linewidth',lw);
ylr = m*T0+q;
h5 = plot(T0,ylr,'g','linewidth',lw);
diffpylr = log(normdu(1))-ylr(1);
qp = log(normdu(1));
h6 = plot(T0,ylr+diffpylr,'r','linewidth',lw);
kk = 750;
diffmylr = log(normdu(kk))-ylr(kk);
h7 = plot(T0,ylr+diffmylr,'b','linewidth',lw);
qm = ylr(1)+diffmylr(1);
output2 = [m exp([q qp qm])]';

h1 = plot(T0,log(normdu),'k','linewidth',lw);
h2 = plot([T0(1) T0(end)], log(avnormdu)*[1 1],'k--','linewidth',lw);

set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'xtick',0:1:T0(end),'fontsize',ftsz)
legend([h1 h6 h5 h7],{'$\ln(RMSE_{u}(t))$ [ms$^{-2}$]',...
    '$\ln(\varepsilon_{u}^{+}(t))$ [ms$^{-2}$]',...
    '$\ln(\varepsilon_{u}(t))$ [ms$^{-2}$]',...
    '$\ln(\varepsilon_{u}^{-}(t))$ [ms$^{-2}$]'},...
    'FontSize',ftsz,'location','northeast',...
        'interpreter','latex')

end


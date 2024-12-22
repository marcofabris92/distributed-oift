function [] = average_u(T0,Ui,h_IN,t_settl,t_settl_2,t_settl_3,...
    ftsz,lw,R,nAg,DIM)

intval = @(i) (DIM*(i-1)+1:DIM*(i-1)+DIM);
log10_05 = @(xx) (log(0.5*xx)/log(10));

R_normalizer = 0;
for i = 1:nAg
    Rii = R(intval(i),intval(i));
    R_normalizer  = R_normalizer + norm(Rii,'fro')^2;
end
R_normalizer = sqrt(R_normalizer);

T = ceil(length(T0)/1);
input_energysq = zeros(nAg,T);
av_input_energysq = zeros(1,T);
tot_input_energysq = zeros(1,T);
for t = 1:T
    for i = 1:nAg
        ui_t = Ui(t,intval(i))';
        input_energysq(i,t) = ui_t'*R(intval(i),intval(i))*ui_t;
    end
    tot_input_energysq(t) = sum(input_energysq(:,t));
    % quadratic Frobenius-weighted mean
    av_input_energysq(t) = tot_input_energysq(t)/R_normalizer;
end

plots = [];
descr = [];
pink = [255,20,147]/255;
puav = plot(T0(1:T),log10_05(av_input_energysq),'g','linewidth',lw);
putot = plot(T0(1:T),log10_05(tot_input_energysq),'color',pink,'linewidth',lw);


av_enL2 = log10_05(trapz(T0(1:T),av_input_energysq)/T0(T));
plot([T0(1) T0(T)],[av_enL2 av_enL2],'g--','linewidth',lw/2)
tot_enL2 = log10_05(2/T0(T)*h_IN);
plot([T0(1) T0(T)],[tot_enL2 tot_enL2],'--','color',pink,'linewidth',lw/2)

plots = [plots puav putot];
descr = [descr, {'$\log_{10}(\bar{l}^{in}(t))$'}, ...
                {'$\log_{10}(l^{in}(t))$'}];

darkgreen = [0 51 25]/255; 
lightgreen = [102 255 178]/255;

% for i = 1:nAg
%     coeff_i = (i-1)/(nAg-1);
%     green_i = lightgreen*(1-coeff_i)+darkgreen*coeff_i;
%     pu = plot(T0,log100(input_energysq(i,:)),'color',green_i,'linewidth',lw/2);
%     av_energy_L2_i = log100(trapz(T0,input_energysq(i,:))/T0(T));
%     plot([T0(1) T0(T)],[av_energy_L2_i av_energy_L2_i],'--','color',green_i,'linewidth',lw/2);
%     plots = [plots pu];
%     str_i = strcat('$E_',num2str(i),'^{in}(t)$');
%     descr = [descr, {str_i}];
% end

max_u = log10_05(max(tot_input_energysq));
max_u = 7;
min1_u = min(min(input_energysq));
min2_u = min(tot_input_energysq);
min3_u = min(av_input_energysq);
min_u = -7; % log10(min([min1_u min2_u min3_u]));

settplot = [];
settplot2 = [];
settplot3 = [];
if t_settl <= T
    settplot = plot([T0(t_settl) T0(t_settl)],[min_u max_u],'k');
end
if t_settl_2 > t_settl && t_settl_2 <= T
    settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[min_u max_u],'--k');
end
if t_settl_3 > t_settl_2 && t_settl_3 <= T
    settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[min_u max_u],'-.k');
end


legend(plots,descr,'FontSize',ftsz,'location','northeast',...
    'interpreter','latex','orientation','horizontal')

xlabel('$t$ [s]','FontSize',ftsz,'interpreter','latex')
ylim([min_u max_u])
%ylabel('$\log_{10}(l^{in}(t))$','FontSize',ftsz,'interpreter','latex')
set(gca,'xtick',0:1:T,'ytick',min_u:1:max_u,'FontSize',ftsz,'TickLabelInterpreter','latex')


end


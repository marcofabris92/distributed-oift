function [] = graphics_main(grp,version,nAg,DIM,ftsz,lw,dijs,R,kr,ka,A)

T0 = grp.T0;
T = length(T0);

XBi_norm = grp.XBi_norm;
XBdes = grp.XBdes;
Ui_norm = grp.Ui_norm;
Xi = grp.Xi;
FOi = grp.FOi;
Ui = grp.Ui;
kT = grp.kT;
gCost_tr_vals = grp.gCost_tr_vals;
gCost_in_vals = grp.gCost_in_vals;
gCost_fo_vals = grp.gCost_fo_vals;
gCost_vals = grp.gCost_vals;
Descent = grp.Descent;

XB = zeros(length(Xi(:,1)), DIM);
for kkk = 1:DIM
  for kk = 1:nAg
    XB(:,kkk) = XB(:,kkk) + Xi(:,(kk-1)*DIM+kkk)/nAg;
  end
end

XBdot = zeros(length(Xi(:,1)), DIM);
for kkk = 1:DIM
  for kk = 1:nAg
    XBdot(:,kkk) = XBdot(:,kkk) + Xi(:,DIM*nAg+(kk-1)*DIM+kkk)/nAg;
  end
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


%% fase portraits
kTgif = 2; % to make animated GIF > 2... if not desired, set it to 2
for ii = 1:kTgif-1
    traj_fig = figure;

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
          txt = '';
          if kt == 1 || kt == kT
            txt = strcat('$',num2str(k+1),'$');
          end
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
          txt = '';
          if kt == 1 || kt == kT
            txt = strcat('$',num2str(k+1),'$');
          end
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
          if DIM == 2 && A(i+1,j+1) == 1
            plot(segment(1,:), segment(2,:), 'color',...
              color_ini/(kT-kt+1)+color_fin*(1-1/(kT-kt+1)))
          elseif DIM == 3 && A(i+1,j+1) == 1
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
if version == 'centralized'
    figure
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
end

%% cost descent (1st order)
if version == 'centralized'
    figure
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
end


% cut the x-axis in a half
% -------------------------------------------------------------------
T = ceil(length(T0)/1);
XBi_norm = grp.XBi_norm(1:T);
XBdes = grp.XBdes(1:T,:);
Xi = grp.Xi(1:T,:);
FOi = grp.FOi(1:T);
Ui = grp.Ui(1:T,:);
XB = XB(1:T,:);
XBdot = XBdot(1:T,:);
% -------------------------------------------------------------------


%% tracking vs formation regimes
figure
grid on
hold on
[t_settl,t_settl_2,t_settl_3,h_TR,h_FO] = tr_vs_fo(T0,XBi_norm,...
    gCost_tr_vals(end),FOi,gCost_fo_vals(end),...
    1,ftsz,lw);
average_tracking_time = -h_TR*T0(end)/(h_FO-h_TR);
average_formation_time = h_FO*T0(end)/(h_FO-h_TR);
fprintf('Average tracking time vs average formation time:\n')
display([average_tracking_time average_formation_time])

%% average & total control for each agent
figure
grid on
hold on
average_u(T0,Ui,gCost_in_vals(end),t_settl,t_settl_2,t_settl_3,...
    ftsz,lw/3*2,R,nAg,DIM);



%% evolutions of the distance errors (sigma_dij-s)
if version == 'distributed'
    M = DIM;
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
                    p_i_t = Xi(tt,(i-1)*M+1:(i-1)*M+M);
                    p_j_t = Xi(tt,(j-1)*M+1:(j-1)*M+M);
                    sij_t = norm(p_i_t-p_j_t)^2;
                    dist_errors(tt,k+1) = sigma(sij_t,dij,0,kr,ka);
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
            if dijs(i,j) > 0
                plot(T0(1:T),dist_errors(:,k+1),tratto,'linewidth',lw/3*2)
            end
            k = k+1;
        end
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],[0 15000],'k'); % ylim
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[0 15000],'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[0 15000],'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','northeast',...
        'interpreter','latex')

    % title('zero-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\sigma(s_{ij}(t))$','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
end


%% evolutions of the first derivatives of the distance errors (sigma_dij-s)

if version == 'distributed'
    M = DIM;
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
                    p_i_t = Xi(tt,(i-1)*M+1:(i-1)*M+M);
                    p_j_t = Xi(tt,(j-1)*M+1:(j-1)*M+M);
                    sij_t = norm(p_i_t-p_j_t)^2;
                    ddist_errors(tt,k+1) = sigma(sij_t,dij,1,kr,ka);
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
            if dijs(i,j) > 0
                plot(T0(1:T),ddist_errors(:,k+1),tratto,'linewidth',lw/3*2)
            end
            k = k+1;
        end
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],[-10 30],'k');
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[-10 30],'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[-10 30],'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','northeast',...
    'interpreter','latex')

    % title('first-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\frac{\partial\sigma}{\partial s_{ij}}(t)$ [m$^{-2}$]','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
end

%% evolutions of the second derivatives of the distance errors (sigma_dij-s)

if version == 'distributed'
    M = DIM;
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
                    p_i_t = Xi(tt,(i-1)*M+1:(i-1)*M+M);
                    p_j_t = Xi(tt,(j-1)*M+1:(j-1)*M+M);
                    sij_t = norm(p_i_t-p_j_t)^2;
                    dddist_errors(tt,k+1) = sigma(sij_t,dij,2,kr,ka);
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
            plot(T0(1:T),dddist_errors(:,k+1),tratto,'linewidth',lw/3*2)
            k = k+1;
        end
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],[0 0.8],'k');
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[0 0.8],'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[0 0.8],'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','northeast',...
        'interpreter','latex')
    
    % title('second-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\frac{\partial^2\sigma}{\partial s_{ij}^2}(t)$ [m$^{-4}$]','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
end


%% evolutions of the relative velocities

if version == 'distributed'
    M = DIM;
    NI = nAg*M;
    figure
    grid on
    hold on

    k = 0;
    for i = 2:nAg
        for j = 1:i-1
            dij = dijs(i,j);
            if dij >= 0
                for tt = 1:T
                    p_dot_i_t = Xi(tt,NI+(i-1)*M+1:NI+(i-1)*M+M);
                    p_dot_j_t = Xi(tt,NI+(j-1)*M+1:NI+(j-1)*M+M);
                    e_dot_ij_t = (p_dot_i_t-p_dot_j_t)';
                    velij_errors(tt,k+1) = norm(e_dot_ij_t);
                end
            else
                velij_errors(:,k+1) = zeros(T,1);
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
            plot(T0(1:T),velij_errors(:,k+1),tratto,'linewidth',lw/3*2)
            k = k+1;
        end
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],ylim,'k');
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],ylim,'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],ylim,'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','northeast',...
        'interpreter','latex')
    
    % title('second-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\left\|\dot{\mathbf{e}}_{ij}(t)\right\|_{2}$ [ms$^{-1}$]','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
end


%% evolutions of the centroid positions

if version == 'distributed'
    figure
    grid on
    hold on

    for k = 1:DIM
        tratto = '-';
        plot(T0(1:T),XB(:,k),tratto,'linewidth',1.5)
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],[-10 30],'k'); % ylim
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[-10 30],'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[-10 30],'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','southeast',...
        'interpreter','latex')
    
    % title('second-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\mathbf{p}_{c}(t)$ [m]','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
    ylim([-10 30])
end


%% evolutions of the centroid velocities

if version == 'distributed'
    figure
    grid on
    hold on

    for k = 1:DIM
        tratto = '-';
        plot(T0(1:T),XBdot(:,k),tratto,'linewidth',lw/3*2)
    end
    
    % settling
    settplot = [];
    settplot2 = [];
    settplot3 = [];
    plots = [];
    descr = [];
    if t_settl <= T
        settplot = plot([T0(t_settl) T0(t_settl)],ylim,'k');
        plots = [plots settplot];
        descr = [descr, {'$10\%$ settling'}];
    end
    if t_settl_2 > t_settl && t_settl_2 <= T
        settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],ylim,'--k');
        plots = [plots settplot2];
        descr = [descr, {'$1\%$ settling'}];
    end
    if t_settl_3 > t_settl_2 && t_settl_3 <= T
        settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],ylim,'-.k');
        plots = [plots settplot3];
        descr = [descr, {'$0.1\%$ settling'}];
    end
    legend(plots,descr,'FontSize',ftsz,'location','southeast',...
        'interpreter','latex')
    
    % title('second-order consensus')
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\dot{\mathbf{p}}_{c}(t)$ [ms$^{-1}$]','interpreter','latex')
    set(gca,'xtick',0:1:T,'fontsize',ftsz)
    %set(get(gca,'ylabel'),'rotation',0)
end

end
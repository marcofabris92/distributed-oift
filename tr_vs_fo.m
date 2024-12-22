function [t_settl,t_settl_2,t_settl_3,h_TR,h_FO] = ...
    tr_vs_fo(T0,l_TR,h_TR,l_FO,h_FO,feas,ftsz,lw)


%fi = @(xx) exp(0.2*log(abs(xx)));


if ~feas
    l_FO = l_FO-min(l_FO);
end
TT = ceil(length(T0(:,1))/1);
max_val = max(abs(l_TR-l_FO));
eps = 10^-50;
values = (l_FO-l_TR)./(1+l_FO+l_TR);
settl_vals = values;
if eps > 10^(-17)
    for t = 1:TT
        if abs(values(t)) < eps || isnan(values(t))
            values(t) = 0;
        end
    end
end
h_TR = -(2/T0(TT))*h_TR/max_val;
h_FO = (2/T0(TT))*h_FO/max_val;
h_TRFO = h_TR+h_FO;

marker_sign = 0;
pos_set = [];
zero_set = [];
neg_set = [];
pos_val = [];
zero_val = [];
neg_val = [];
red_trigger = 0;
blue_trigger = 0;
yellow_trigger = 0;
red_plot = [];
blue_plot = [];
yellow_plot = [];
for t = 1:TT
    new_sign = sign(values(t));
    switch new_sign 
        case 1
            pos_set = [pos_set T0(t)];
            pos_val = [pos_val values(t)];
        case -1
            neg_set = [neg_set T0(t)];
            neg_val = [neg_val values(t)];
        otherwise 
            zero_set = [zero_set T0(t)];
            zero_val = [zero_val values(t)];
    end
    if t == TT || t > 1 && marker_sign ~= new_sign
        switch marker_sign
            case 1
                if length(pos_set) == 1
                    red_plot = plot([0 0],[-5 -5],'r','linewidth',lw);
                    plot(pos_set,pos_val,'ro','MarkerFaceColor','r','linewidth',lw)
                else
                    red_plot = plot(pos_set,pos_val,'r','linewidth',lw);
                    area(pos_set,pos_val,'FaceColor','r')
                end
                pos_set = [];
                pos_val = [];
            case -1
                if length(neg_set) == 1
                    blue_plot = plot([0 0],[-5 -5],'b','linewidth',lw);
                    plot(neg_set,neg_val,'bo','MarkerFaceColor','b','linewidth',lw)
                else
                    blue_plot = plot(neg_set,neg_val,'b','linewidth',lw);
                    area(neg_set,neg_val,'FaceColor','b')
                end
                neg_set = [];
                neg_val = [];
            otherwise
                
                if length(zero_set) == 1
                    yellow_plot = plot([0 0],[-5 -5],'y','linewidth',lw);
                    plot(zero_set,zero_val,'yo','MarkerFaceColor','y','linewidth',lw)
                else
                    yellow_plot = plot(zero_set,zero_val,'y','linewidth',lw);
                end
                zero_set = [];
                zero_val = [];
        end
    end
    marker_sign = new_sign;
end

plots = [];
descr = {};
if ~isempty(red_plot) 
    plots = [plots red_plot];
    descr = {'Formation regime'};
end
if ~isempty(yellow_plot) 
    plots = [plots yellow_plot];
    if ~isempty(red_plot) 
        descr = {'Formation regime','Balanced regime'};
    end
end
if ~isempty(blue_plot) 
    plots = [plots blue_plot];
    if ~isempty(red_plot)
        if ~isempty(yellow_plot)
            descr = {'Formation regime','Balanced regime','Tracking regime'};
        else
            descr = {'Formation regime','Tracking regime'};
        end
    else
        if ~isempty(yellow_plot)
            descr = {'Balanced regime','Tracking regime'};
        else
            descr = {'Tracking regime'};
        end
    end
end

yellow = [1 1 0];
ltrfoborder = plot(T0(1:TT),values,'color',yellow,'linewidth',lw/2);

% average plots
%plot([0 T0(TT)],[h_TR h_TR],'--b','linewidth',lw/2)
%plot([0 T0(TT)],[h_FO h_FO],'--r','linewidth',lw/2)
%plot([0 T0(TT)],[h_TRFO h_TRFO],'--','color',yellow,'linewidth',lw/2);



%% settling plots
settling = 0.1;
plot([0 T0(TT)],[settling settling],'k')
plot([0 T0(TT)],[-settling -settling],'k')
plot([0 T0(TT)],[settling^2 settling^2],'--k')
plot([0 T0(TT)],[-settling^2 -settling^2],'--k')
%plot([0 T0(TT)],[settling^3 settling^2],'--k')
%plot([0 T0(TT)],[-settling^3 -settling^2],'--k')

t = TT;
trigger_settl_2 = 0;
trigger_settl_3 = 0;
done = 0;
t_settl = TT+1;
t_settl_2 = TT+1;
t_settl_3 = TT+1;
while ~done && t > 0
    if abs(settl_vals(t)) > settling^3 && ~trigger_settl_3
        trigger_settl_3 = 1;
    elseif ~trigger_settl_3
        t_settl_3 = t;
    end
    if abs(settl_vals(t)) > settling^2 && ~trigger_settl_2
        trigger_settl_2 = 1;
    elseif ~trigger_settl_2
        t_settl_2 = t;
    end
    if abs(settl_vals(t)) > settling
        done = 1;
    else
        t_settl = t;
    end
    t = t-1;
end
settplot = [];
settplot2 = [];
settplot3 = [];
plots = [plots ltrfoborder];
descr = [descr, {'$l^{tf}(t)$'}];
if t_settl <= TT
    settplot = plot([T0(t_settl) T0(t_settl)],[-1 1],'k');
    plots = [plots settplot];
    descr = [descr, {'$10\%$ settling'}];
end
if t_settl_2 > t_settl && t_settl_2 <= TT
    settplot2 = plot([T0(t_settl_2) T0(t_settl_2)],[-1 1],'--k');
    plots = [plots settplot2];
    descr = [descr, {'$1\%$ settling'}];
end
if t_settl_3 > t_settl_2 && t_settl_3 <= TT
    settplot3 = plot([T0(t_settl_3) T0(t_settl_3)],[-1 1],'-.k');
    plots = [plots settplot3];
    descr = [descr, {'$0.1\%$ settling'}];
end



legend(plots,descr,'FontSize',27,'location','northeast','interpreter','latex')

xlabel('$t$ [s]','FontSize',ftsz,'interpreter','latex')
ylim([-1 1])
ylabel('$l^{tf}(t)$','FontSize',ftsz,'interpreter','latex')
set(gca,'xtick',0:1:TT,'ytick',-1:0.2:1,'FontSize',ftsz,'TickLabelInterpreter','latex')

end
function cG = graphics_est(G,NI,M,kpy,kdy,xx_history,u_history,T0,ftsz,lw)

NS = 2*NI;
K = length(T0);
n = numnodes(G);
L = full(laplacian(G));
NN = zeros(NI*n);
IM = eye(M);
for i = 1:n
    index_i = (i-1)*NI+(i-1)*M+1:(i-1)*NI+(i-1)*M+M;
    NN(index_i,index_i) = IM;
end
MM = (kron(L,eye(NI))+NN)^-1;
cG = max(eig(MM));

NIE = NI*n;
p_history = xx_history(:,1:NI);
pdot_history = xx_history(:,NI+1:NS);
n2epc = zeros(K,1);
n2epdotc = zeros(K,1);
UBp = zeros(K,1);
UBpdot = zeros(K,1);
for k = 1:K
    for i = 1:n
        
        epcki = zeros(1,M);
        epdotcki = zeros(1,M);
        for j = 1:n
            index_j = (j-1)*M+1:(j-1)*M+M;
            index_pij = NS+(i-1)*NI+index_j;
            index_pdotij = NS+NIE+(i-1)*NI+index_j;
            epcki = epcki + p_history(k,index_j)-xx_history(k,index_pij);
            epdotcki = epdotcki + pdot_history(k,index_j)-...
                xx_history(k,index_pdotij);
        end
        epcki = epcki/n;
        epdotcki = epdotcki/n;
        n2epc(k) = n2epc(k) + norm(epcki)^2;
        n2epdotc(k) = n2epdotc(k) + norm(epdotcki)^2;
    
        Ni = neighbors(G,i);
        Nicompl = setdiff((1:n),union(Ni,i));
        wipk = 0;
        wipdotk = 0;
        for j_ = 1:length(Nicompl)
            j = Nicompl(j_);
            
            index_wjk = (j-1)*M+1:(j-1)*M+M;
            index_wtilijk = NS+NIE+(i-1)*NI+index_wjk;
            dwipk = pdot_history(k,index_wjk)-...
                xx_history(k,index_wtilijk);
            dwipdotk = u_history(k,index_wjk);
            wipk = wipk + norm(dwipk)^2;
            wipdotk = wipdotk + norm(dwipdotk)^2;
        end
        
        UBp(k) = UBp(k) + wipk;
        UBpdot(k) = UBpdot(k) + wipdotk;
    end
end
UBp = UBp*(cG/(n*kpy))^2;
UBpdot = UBpdot*(cG/(n*kdy))^2;


figure
grid on
hold on
h1 = plot(T0,log10(n2epc),'m','linewidth',lw);
h2 = plot(T0,log10(n2epdotc),'c','linewidth',lw);
h3 = plot (T0,log10(UBp),'r','linewidth',lw);
h4 = plot (T0,log10(UBpdot),'b','linewidth',lw);
% plot(T0,log10(n2epc+n2epdotc),'r')
% plot (T0,log10(UBp+UBpdot),'b')

set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'xtick',0:1:T0(end),'fontsize',ftsz)
legend([h1 h2 h3 h4],{'$\log_{10}(\left\| \mathbf{e}_{p_{c}}(t) \right\|_{2}^{2})$ [m]',...
    '$\log_{10}(\left\| \mathbf{e}_{\dot{p}_{c}}(t) \right\|_{2}^{2})$ [ms$^{-1}$]',...
    '$\log_{10}(\varpi_{p_{c}}^{2}(t))$ [m]',...
    '$\log_{10}(\varpi_{\dot{p}_{c}}^{2}(t))$ [ms$^{-1}$]'},...
    'FontSize',ftsz,'location','northeast',...
        'interpreter','latex')

end


function [gCost_tri, gCost_ini, gCost_foi, XBi_norm, Ui_norm, FOi] =...
    computeTrInFoCosts(Xi,Ui,Xdes,QB,R,nAg,DIM,T0,dijs,kF,kA,kr,ka)


NI = nAg*DIM;

TT = length(T0);
XBi_norm = zeros(TT,1);   % nAg/2*|| xB - xB,des ||_QB ^2 for all t in [0,T]
Ui_norm = zeros(TT,1);    % 1/2*|| u ||_R ^2 for all t in [0,T]
FOi = zeros(TT,1);        % kF/2*(sum_1^n sum_j\in Ni sigma(sij))
DxBi_t = zeros(2*DIM,1);  % xB-xB,des at time t

for t = 1:TT
    for j = 1:2*DIM
        DxBi_t(j) = 0;
    end
    for j = 1:DIM
        for jj = 1:nAg
            k = (jj-1)*DIM+j;
            DxBi_t(j) = DxBi_t(j) + (Xi(t,k) - Xdes(t,j));
            DxBi_t(DIM+j) = DxBi_t(DIM+j) + (Xi(t,NI+k) - Xdes(t,DIM+j));
            
        end
    end
    for j = 1:2*DIM
        DxBi_t(j) = DxBi_t(j)/nAg;
    end
    XBi_norm(t) = nAg*DxBi_t'*QB*DxBi_t/2;
    Ui_norm(t) = Ui(t,1:NI)*R*Ui(t,1:NI)'/2;
    FOi(t) = l_FO(Xi(t,1:NI),Xi(t,NI+1:2*NI),nAg,DIM,dijs,kF,kA,kr,ka);
end

gCost_tri = trapz(T0,XBi_norm);
gCost_ini = trapz(T0,Ui_norm);
gCost_foi = trapz(T0,FOi);
  
end


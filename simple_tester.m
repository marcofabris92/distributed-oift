clear all, clc

x_des = [0.0, 0.0, 0.0;   0.0, 0.0, 0.0;   0.0, 0.0, 0.0; ...
        1.0, 0.0, 0.0;   1.0, 0.0, 0.0;   1.0, 0.0, 0.0];
% x = [0 0.5 0; 0 -0.5 0; sqrt(3)/2 0 0; 1 0 0; 1 0 0; 1 0 0];
x = [0.0, 0.0, 0.0;  -2.0, 1.0, 0.0;   0.0, -1.0, 0.0;...
    1.0, 0.0, 0.0;   1.0, 0.0, 0.0;   1.0, 0.0, 0.0];
u = [0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0 ];
 
nAg = 3;
DIM = 3;
dd = 1;
I = eye(DIM);
H_p_p_FF = zeros(nAg*DIM);
nabla_p_FF = zeros(nAg,DIM);
do_Qsafe = 0;
scale = 1/( dd*(sqrt(2)-1) )^2;
for k = 1:nAg
    for l = 1:nAg
        
        if k ~= l

            p_i = x(k,:)';
            p_j = x(l,:)';
            del_p_ij = p_i-p_j;
            sij = norm(del_p_ij)^2;
            Pij = del_p_ij*del_p_ij';

            Dsig_ij = scale*sign(sij-dd)*( 1-dd/(abs(sij-dd^2)+dd^2)^(1/2) );
            DDsig_ij = scale*dd/(2*(abs(sij-dd^2)+dd^2)^(3/2));
            
            kki = (k-1)*DIM+1;
            kkf = kki+DIM-1;
            lli = (l-1)*DIM+1;
            llf = lli+DIM-1;
            
            H_p_p_FF(kki:kkf,lli:llf) = -4*DDsig_ij*Pij;
            if Dsig_ij >= 0
                H_p_p_FF(kki:kkf,lli:llf) = H_p_p_FF(kki:kkf,lli:llf)-2*Dsig_ij*I;
            else
                do_Qsafe = 1;
            end
            H_p_p_FF(kki:kkf,kki:kkf) = H_p_p_FF(kki:kkf,kki:kkf)-H_p_p_FF(kki:kkf,lli:llf);
            
            nabla_p_FF(k,:) = nabla_p_FF(k,:) + Dsig_ij*(2*del_p_ij)';
        end
            
    end
end

nabla_p_FF
H_p_p_FF
eigH = eig(H_p_p_FF)
if do_Qsafe
    fprintf('DO Qsafe!\n')
end


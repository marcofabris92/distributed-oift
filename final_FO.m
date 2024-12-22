function result = final_FO(Xf,n,DIM,nAg,kF,kA,dijs,kr,ka)

NI = DIM*nAg;
NS = 2*NI;
I_DIM = eye(DIM);

p = zeros(DIM,nAg);
dp = zeros(DIM,nAg);
for i = 1:nAg
    p(:,i) = Xf((i-1)*DIM+1:(i-1)*DIM+DIM)';
    dp(:,i) = Xf(NI+(i-1)*DIM+1:NI+(i-1)*DIM+DIM)';
end

if n == 0
    cost_FO = 0;
    for i = 1:nAg
        for j = 1:i-1
            eij = p(:,i)-p(:,j);
            deij = dp(:,i)-dp(:,j);
            sij = eij'*eij;
        	dij = dijs(i,j);
            if dij >= 0
                cost_FO = cost_FO + kF*sigma(sij,dij,0,kr,ka) + kA*(deij'*deij);
            end
        end
    end
    result = cost_FO/2;
end

if n == 1
    grad_FO = zeros(NS,1);
    for i = 1:nAg
        intval_i = (i-1)*DIM+1:(i-1)*DIM+DIM;
        for j = 1:nAg
            if j ~= i
                eij = p(:,i)-p(:,j);
                deij = dp(:,i)-dp(:,j);
                sij = eij'*eij;
                dij = dijs(i,j);
                if dij >= 0
                    grad_FO(intval_i) = grad_FO(intval_i) + ...
                        kF*sigma(sij,dij,1,kr,ka)*eij;
                    grad_FO(NI+intval_i) = grad_FO(NI+intval_i) + kA*deij;
                end
            end
        end
    end
    result = grad_FO;
end

if n == 2
    H_FO = zeros(NS);
    for i = 1:nAg
        iint = (i-1)*DIM+1:(i-1)*DIM+DIM;
        for j = 1:nAg
            if j~=i
                dij = dijs(i,j); 
                if dij >= 0
                    jint = (j-1)*DIM+1:(j-1)*DIM+DIM;
                    eij = p(:,i)-p(:,j);
                    PIij = eij*eij';
                    sij = eij'*eij;
                    H_FO(iint,jint) = -kF*2*sigma(sij,dij,2,kr,ka)*PIij;
                    H_FO(NI+iint,NI+jint) = -kA*I_DIM;
                    sig1ij = sigma(sij,dij,1,kr,ka);
                    if sig1ij > 0
                        H_FO(iint,jint) = H_FO(iint,jint) - kF*sig1ij*I_DIM;
                    end
                    H_FO(iint,iint) = H_FO(iint,iint) - H_FO(iint,jint);
                    H_FO(NI+iint,NI+iint) = H_FO(NI+iint,NI+iint)...
                        - H_FO(NI+iint,NI+jint);
                end
            end
        end
    end
    
    % eigenvalue compensation to get H_FO pos. semi-def.
%     specH_FO = sort(eig(H_FO(1:NS/2,1:NS/2)));
%     lambdaH_FO_0 = abs(specH_FO(1));
%     for k = 1:NS/2
%         H_FO(k,k) = H_FO(k,k) + lambdaH_FO_0;
%     end
    
    result = H_FO;
end



% result = 1*result;

end
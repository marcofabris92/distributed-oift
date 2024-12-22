function [] = rigidity_analysis(G,p,nAg,M)


Nedges = numedges(G);
diageij = zeros(Nedges,Nedges*M);
E = zeros(nAg,Nedges);
k = 0;
for i = 2:nAg
    for j = 1:i-1
        diageij(k+1,k*M+1:k*M+M) = p((j-1)*M+1:(j-1)*M+M)-...
            p((i-1)*M+1:(i-1)*M+M);
        E(i,k+1) = -1;
        E(j,k+1) = 1;
        k = k+1;
    end
end
Rig = diageij*kron(E',eye(M));
rkRig = rank(Rig);
target_rank = length(p)-(M*(M+1))/2;

% diageij
% E
% Rig
% target_rank
% rkRig
% Nedges

fprintf(' ')
fprintf(num2str(target_rank))
fprintf(' ')
fprintf(num2str(rkRig))
fprintf(' ')
fprintf(num2str(Nedges))
fprintf(' ')
if target_rank == rkRig && target_rank == Nedges
    fprintf('IMR\n')
elseif target_rank == rkRig
    fprintf('IR\n')
end


end


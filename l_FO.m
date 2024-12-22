function value = l_FO(p,dp,nAg,DIM,dijs,kF,kA,kr,ka)

if kF <= 0
    error('\nkF must be positive\n')
end

value = 0;

for i = 2:nAg
    for j = 1:i-1
        if dijs(i,j) >= 0
            sij = norm(p((i-1)*DIM+1:(i-1)*DIM+DIM)-...
                p((j-1)*DIM+1:(j-1)*DIM+DIM))^2;
            alij = norm(dp((i-1)*DIM+1:(i-1)*DIM+DIM)-...
                dp((j-1)*DIM+1:(j-1)*DIM+DIM))^2;
            value = value + kF*sigma(sij,dijs(i,j),0,kr,ka) + kA*alij;
        end
    end
end

value = value/2;

end
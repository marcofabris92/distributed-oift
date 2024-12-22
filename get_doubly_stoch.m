% metropolis-hastings rule to get a doubly stoch. matrix

function P = get_doubly_stoch(A)

L = length(A);
P = A;
degs = sum((A'));
for i = 1:L
    pijs = 0;
    for j = 1:L
        if i ~= j
            P(i,j) = 1/(1+max(degs(i),degs(j)));
            pijs = pijs + P(i,j);
        end
    end
    P(i,i) = 1-pijs;
end

end


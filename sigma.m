function value = sigma(sij,dij,n,kr,ka)

if sij < 0 || dij < 0 || kr <= 0 || ka <= 0
    error('\nInput quantities must be positive or nonnegative.\n')
end

switch n
    case 0
        if sij <= dij^2
            value = kr*(1-sij/dij^2)^3;
        else
            value = ka*(sqrt(sij)/dij-1)^3;
        end
    case 1
        if sij <= dij^2
            value = 3*kr*(1-sij/dij^2)^2*(-1/dij^2);
        else
            value = 3*ka*(sqrt(sij)/dij-1)^2*(1/(2*sqrt(sij)*dij));
        end
    case 2
        if sij <= dij^2
            value = 6*kr*(1-sij/dij^2)*(1/dij^4);
        else
            value = (3*ka*(sij-dij^2))/(4*dij^4*sij^1.5);
        end
end

end
function [a_hat,E] = ARmodelR1a(x,N)
    y = x(N+1:end);
    H = x(1:end-N);
    a_hat = (H'*H)\H'*y;
    E = norm(y-H*a_hat)^2;
end


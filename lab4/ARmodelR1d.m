function [a_hat,E] = ARmodelR1d(r,P)
    y = r(P+1:end);
    H = zeros(length(r)-P,P);
    for i = P:length(r)-1
       H(i-P+1,:) = r(i:-1:i+1-P)'; 
    end
    a_hat = (H'*H)\H'*y;
    E = norm(y-H*a_hat)^2;
end


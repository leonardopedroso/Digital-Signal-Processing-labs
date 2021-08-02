function [mu,var] = MLEnormal(x)
    mu = sum(x)/length(x);
    var = sum((x-mu).^2)/length(x);
end
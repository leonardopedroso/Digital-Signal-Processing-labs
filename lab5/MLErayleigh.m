function f = MLErayleigh(x)
    f = sqrt(sum(x.^2)/(2*length(x)));
end
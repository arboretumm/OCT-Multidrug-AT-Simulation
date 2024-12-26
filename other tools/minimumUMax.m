function [minUMax] = minimumUMax(p)
    gamma = p.lambda_n;
    kappaR = p.kappa_threshold/p.kappa;
    alphaA = p.alphaA;
    minUMax = (gamma.*(1-kappaR))./alphaA;
end
function [worstTime] = calculateWorstFromPopStruct(p)
    kappa_r = p.kappa_threshold/p.kappa;
    kappa_d = p.doubleMutPopIC/p.kappa;
    worstTime = (-1/p.lambda_d).*log(((p.kappa - p.kappa_threshold)/p.kappa_threshold)*(1/((1/kappa_d)-1)));
end
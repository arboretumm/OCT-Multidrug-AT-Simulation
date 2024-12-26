function [idealTime] = calculateIdealFromPopStruct(p)
    kappa_r = p.kappa_threshold/p.kappa;
    kappa_d = p.doubleMutPopIC/p.kappa;
    idealTime = (1./(p.lambda_d.*(1-kappa_r)).*log(kappa_r/kappa_d));
end
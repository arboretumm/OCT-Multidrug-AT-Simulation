function [position, isterminal, direction] = eventMaxB_Anaught(t, x, p)
    concB = ((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* (1 - (p.kappa_threshold)/p.kappa))/(p.alphaB .* (x(1) + x(3)));
    position = concB - p.uMaxB;
    isterminal = 1;
    direction = 1;
end
function [position, isterminal, direction] = eventMaxA_Bmax(t, x, p)
    concA = (((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* (1 - (p.kappa_threshold)/p.kappa)) - p.uMaxB.*p.alphaB.*(x(1) + x(3)))/(p.alphaA .* (x(1) + x(2)));
    position = concA - p.uMaxA;
    isterminal = 1;
    direction = 1;
end
function [position, isterminal, direction] = eventGoalA_BSub(t, x, p, goal, other)
    concA = (((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* (1 - (p.kappa_threshold)/p.kappa)) - other.*p.alphaB.*(x(1) + x(3)))/(p.alphaA .* (x(1) + x(2)));
    position = concA - goal;
    isterminal = 1;
    direction = 1;
end
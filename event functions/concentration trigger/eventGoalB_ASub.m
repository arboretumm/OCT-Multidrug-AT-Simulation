function [position, isterminal, direction] = eventGoalB_ASub(t, x, p, goal, other)
    %only for threshold concentration
    concB = ((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* (1 - (p.kappa_threshold)/p.kappa) - p.alphaA.*other.*(x(1) + x(2)))/(p.alphaB .* (x(1) + x(3)));
    position = concB - goal;
    isterminal = 1;
    direction = 1;
end
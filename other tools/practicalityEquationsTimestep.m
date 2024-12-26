function dx = practicalityEquationsTimestep(t, x, p, concA, concB)

dx = [0; 0; 0; 0; 0];

%equations
    dx(1) = p.lambda_n .* x(1) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - (p.alphaA * x(1) * concA) - (p.alphaB * x(1) * concB);
    dx(2) = p.lambda_m .* x(2) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaA * x(2) * concA;
    dx(3) = p.lambda_s .* x(3) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaB * x(3) * concB;
    dx(4) = p.lambda_d .* x(4) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa);
    dx(5) = dx(1) + dx(2) + dx(3) + dx(4);
end
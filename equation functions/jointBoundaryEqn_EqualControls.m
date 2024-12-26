function dx = jointBoundaryEqn_EqualControls(t, x, p)

dx = [0; 0; 0; 0; 0];

concBoth = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)))./(p.alphaA.*(x(1)+x(2))+p.alphaB.*(x(1)+x(3)));

if (x(1)==0 && x(3)==0 && x(2)==0) || concBoth < 0
    concA = 0;
    concB = 0;
elseif concBoth > p.uMaxB && concBoth <= p.uMaxA
    concA = concBoth;
    concB = p.uMaxB;
elseif concBoth > p.uMaxA && concBoth <= p.uMaxB
    concA = p.uMaxA;
    concB = concBoth;
elseif concBoth > p.uMaxA && concBoth > p.uMaxB
    concA = p.uMaxA;
    concB = p.uMaxB;
else
    concA = concBoth;
    concB = concBoth;
end

dx(1) = p.lambda_n .* x(1) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - (p.alphaA * x(1) * concA) - (p.alphaB * x(1) * concB);
dx(2) = p.lambda_m .* x(2) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaA * x(2) * concA;
dx(3) = p.lambda_s .* x(3) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaB * x(3) * concB;
dx(4) = p.lambda_d .* x(4) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa);
dx(5) = dx(1) + dx(2) + dx(3) + dx(4);

end
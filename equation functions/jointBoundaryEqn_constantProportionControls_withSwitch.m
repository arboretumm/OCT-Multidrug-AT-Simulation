function dx = jointBoundaryEqn_constantProportionControls_withSwitch(t, x, p, proportionConst)

dx = [0; 0; 0; 0; 0];

concCalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)))./(p.alphaA.*(x(1)+x(2))+p.alphaB.*proportionConst.*(x(1)+x(3)));
concCalcB = proportionConst.*concCalcA;

%if concCalcA hits min or max, will recalculate B separately
if (x(1)==0 && x(3)==0 && x(2)==0) 
    concA = 0;
    concB = 0;
elseif concCalcA < 0
    concA = 0;
    concCalcB = (1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) - (p.alphaA.*concA.*(x(1)+x(2)))./(p.alphaB.*(x(1)+x(3))); %insertequationusingconcAinstead
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
elseif concCalcA > p.uMaxA
    concA = p.uMaxA;
    concCalcB = (1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) - (p.alphaA.*concA.*(x(1)+x(2)))./(p.alphaB.*(x(1)+x(3))); %insertequationusingconcAinstead
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
else
    concA = concCalcA;
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
end

dx(1) = p.lambda_n .* x(1) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - (p.alphaA * x(1) * concA) - (p.alphaB * x(1) * concB);
dx(2) = p.lambda_m .* x(2) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaA * x(2) * concA;
dx(3) = p.lambda_s .* x(3) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaB * x(3) * concB;
dx(4) = p.lambda_d .* x(4) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa);
dx(5) = dx(1) + dx(2) + dx(3) + dx(4);

end
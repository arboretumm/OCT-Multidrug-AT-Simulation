function dx = elsasEquations(t, x, p)

dx = [0; 0; 0; 0; 0;];

%[gamma*k_T*(1- k_T/k) ] / [ 2*alpha_B*(x_1+x_2) ].
concCalcB = (p.lambda_n.*p.kappa_threshold.*(1 - p.kappa_threshold./p.kappa))./(2.*p.alphaB.*(x(1)+x(2)));
concCalcA = (((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* ((1-(p.kappa_threshold./p.kappa)))) - p.alphaB.*concCalcB.*(x(1) + x(3)))/(p.alphaA .* (x(1) + x(2)));

if (x(1)==0 && x(2)==0) 
    concA = 0;
    concB = 0;
elseif concCalcA < 0
    concA = 0;
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
elseif concCalcA > p.uMaxA
    concA = p.uMaxA;
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
else
    if concCalcB < 0
        concB = 0;
    elseif concCalcB > p.uMaxB
        concB = p.uMaxB;
    else
        concB = concCalcB;
    end
    reconcA = (((p.lambda_n .* x(1) + p.lambda_m .* x(2) + p.lambda_s .* x(3) + p.lambda_d .* x(4)) .* ((1-(p.kappa_threshold./p.kappa)))) - p.alphaB.*concB.*(x(1) + x(3)))/(p.alphaA .* (x(1) + x(2)));

    if reconcA < 0
        concA = 0;
    elseif reconcA > p.uMaxA
        concA = p.uMaxA;
    else
        concA = reconcA;
    end
end

dx(1) = p.lambda_n .* x(1) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - (p.alphaA * x(1) * concA) - (p.alphaB * x(1) * concB);
dx(2) = p.lambda_m .* x(2) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaA * x(2) * concA;
dx(3) = p.lambda_s .* x(3) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa) - p.alphaB * x(3) * concB;
dx(4) = p.lambda_d .* x(4) .* (1 - (x(1)+x(2)+x(3)+x(4))/p.kappa);
dx(5) = dx(1) + dx(2) + dx(3) + dx(4);

end
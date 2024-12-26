%variables i want to play with
    %alphaVal = 0.00155; %10nM val
    alphaVal = 0.0005175; %30nM val
    ic50Conc = 30;

%variables that are set
    initCondX1 = 40000;
    gammaVal = 0.031;
    kappa = initCondX1;


tSpan = linspace(0, 1000, 100);
[t, x] = ode45(@(t, x) gammaVal*x*(1-x/kappa) - (alphaVal*x*ic50Conc), tSpan, initCondX1);

figure(1)
hold on;
plot(t, x, 'LineWidth', 3)
yline((0.5*initCondX1), 'LineWidth', 2, 'LineStyle','--', 'Color', '#ccc')
xlabel("hours", 'FontSize', 18)
ylabel("cell number", 'FontSize', 18)
axis([0 inf 0 inf])
title(sprintf("Alpha Val to IC50 comp - %s for alpha", num2str(alphaVal)));
set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 16)
hold off;


% function dx = justNaive(t, x)
%     dx = 0;
%     dx(1) = gammaVal* x(1) .* (1 - x(1)/kappa) - (alphaVal * x(1) * ic50Conc);
% end
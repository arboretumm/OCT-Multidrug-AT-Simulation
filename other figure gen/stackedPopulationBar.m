x = [1, 2, 3, 4, 5];
y = zeros(5, 4);
ratio_2naught_3naught = [.05, .25, .5, .75, .95];

for i=1:5
    y(i, 1) = 68.5;
    y(i, 3) = ratio_2naught_3naught(i)*30;
    y(i, 2) = (1 - ratio_2naught_3naught(i))*30;
    y(i, 4) = 1.5;
end

figure(1)
    bar(x, y, 'stacked', 'BarWidth', .6)
    xlabel('different population initial conditions')
    ylabel('% of initial population')
    legend({'x_1', 'x_2', 'x_3', 'x_4'})
    set(gcf, 'units', 'points', 'position', [10, 10, 810, 648]);

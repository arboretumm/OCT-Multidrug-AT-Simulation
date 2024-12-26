function [position, isterminal, direction] = eventSubpopulationLimit(t, x, p, populationNum, populationLim)
    %triggers a stop if a subpopulation hits a particular threshold
    position = x(populationNum) - populationLim;
    isterminal = 1;
    direction = 1;
end
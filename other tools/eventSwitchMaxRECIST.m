function [position, isterminal, direction] = eventSwitchMaxRECIST(t, x, p, xMin)
    position = x(5) - 1.2.*xMin;
    isterminal = 1;
    direction = 1;
end
function [position, isterminal, direction] = eventFindTreatmentMinima_PracticalityA(t, x, p, timeVector, concAVec, concBVec)
    dxdt = practicalityEquations(t, x, p, timeVector, concAVec, concBVec);
    position = dxdt(5);
    isterminal = 1;
    direction = 1;
end
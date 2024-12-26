function uMaxA_uSingularB(folderName, fileName, pStruct)
    %population + time parameters
    p = pStruct;
    initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];

    %solving/graphing system given parameters
    tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
    %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);

    Store.t=[];
    Store.x=[];
    t=tSPAN(1);
    x=initCond;

    if t(end) < tSPAN(end)
        options = odeset('Events', @(t,x) eventSubpopulationLimit(t, x, p, 5, (p.kappa_threshold+1)),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
        [t, x] = ode45(@(t, x) controlAMax_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];

        if t(end) < tSPAN(end)
            options3 = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
            [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options3);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
        end
    end

    %GRAPHING
gcf1 = figure(1);
plot(Store.t, Store.x, 'LineWidth', 2.5);
axis([0 tSPAN(end) 0 inf]);
title('Plot of cell populations', 'FontSize', 24);
xlabel('time (hours)', 'FontSize', 18);
ylabel('cell number', 'FontSize', 18);
legend('susceptible', 'B resist', 'A resist', 'both resist', 'total');
set(gcf1, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 16)

strFileName1 = sprintf('%s%s.fig', folderName, fileName);
strFileName2 = sprintf('%s%s.png', folderName, fileName);
saveas(figure(1),strFileName1); %saves population graph as fig
saveas(figure(1),strFileName2); %saves population graph as png

% %parameter file save
% strFileName3 = sprintf('%s%s.mat', folderName, fileName);
% save(strFileName3, 'p', 'c', "runType");
[rowsz, colsz] = size(Store.x);
foundBreachC = false;
controlFailureTrack = 0;

for q=1:rowsz
    if Store.x(q, 4) > p.kappa_threshold && foundBreachC == false
        disp('simulation broke threshold at');
        disp(Store.t(q));
        controlFailureTrack = Store.t(q);
        foundBreachC = true;
    end
end

%parameter file save with post calc
strFileName3 = sprintf('%s%s.mat', folderName, fileName);
save(strFileName3, 'p', 'runType', 'controlFailureTrack');

toc

end
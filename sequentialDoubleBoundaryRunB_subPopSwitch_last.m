function [controlFailureTrack1, controlFailureTrack2] = sequentialDoubleBoundaryRunB_subPopSwitch_last(folderName, fileName, pStruct, populationNum, populationLim)
%file and folderName save the [ANYTHING ELSE?] and the figures generated in
%the filepath created 
close all;
%disp('seqB last')
%GENERAL RULES
    %use camelcase, not underscores for parameter names EXCEPTING things
    %that would have otherwise been sub or superscripted or need to split
    %apart chunks of a name
    %dear god comment you blessed idiot

    runType = 'sequentialDoubleBoundaryControlShapeB_subPopSwitch';
    %sequential double bound with A as the first boundary, B as 0 -> second
    %boundary, with switches between different phases based on sub
    %population for the first one, and currently just on when the second
    %bound hits max concentration

    %parameter definition
    %base equation build + instantiation of structure
        %structures now built outside of main run file (this one) and fed
        %in as full structures to this file in order to set up multi sim
        %runs easier

    eventBreakTimes = zeros(1, 2);

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
        %phase 1 - control B at 0 and control A at threshold until
        %time-based switch
        options = odeset('Events', @(t,x) eventSubpopulationLimit(t, x, p, populationNum, populationLim),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
        %sprintf("phase 1 start")
        [t, x] = ode45(@(t, x) controlAZero_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        eventBreakTimes(1) = Store.t(end);
        %sprintf("phase 1 done")
        if t(end) < tSPAN(end)
            %phase 2 - control B at threshold and control A at max until
            %time-based switch
            %sprintf("phase 2 start")
            options2 = odeset('Events', @(t,x) eventSubpopulationLimit(t, x, p, populationNum, populationLim),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
            [t, x] = ode45(@(t, x) controlACalcThreshold_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options2);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
            %sprintf("phase 2 done")
            eventBreakTimes(2) = Store.t(end);
            if t(end) < tSPAN(end)
                %phase 3 - control B and control A at max until end of
                %simulation
                %sprintf("phase 3 start")
                options3 = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options3);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
                %sprintf("phase 3 done")
            end
        end
    end

%GRAPHING
gcf1 = figure('Visible', 'off');
plot(Store.t, Store.x, 'LineWidth', 2.5);
axis([0 tSPAN(end) 0 inf]);
title('Plot of cell populations', 'FontSize', 24);
xlabel('time (hours)', 'FontSize', 18);
ylabel('cell number', 'FontSize', 18);
legend('susceptible', 'B resist', 'A resist', 'both resist', 'total');
set(gcf1, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 16)
set(gcf1, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');

strFileName1 = sprintf('%s%s.fig', folderName, fileName);
strFileName2 = sprintf('%s%s.png', folderName, fileName);
savefig(gcf1,strFileName1); %saves population graph as fig
%saveas(gcf1,strFileName2); %saves population graph as png

% %parameter file save
% strFileName3 = sprintf('%s%s.mat', folderName, fileName);
% save(strFileName3, 'p', 'c', "runType");
[rowsz, colsz] = size(Store.x);
controlFailureTrack1 = 0;
controlFailureTrack2 = 0;
errorBreak = 1.001*(p.kappa_threshold);

cf1Index = find((Store.x(:, 4) <= errorBreak), 1, 'last');
cf2Index = find((Store.x(:, 5) <= errorBreak), 1, 'last');

controlFailureTrack1 = Store.t(cf1Index);
controlFailureTrack2 = Store.t(cf2Index);

%parameter file save with post calc
strFileName3 = sprintf('%s%s.mat', folderName, fileName);
save(strFileName3, 'p', 'runType', 'controlFailureTrack1', 'controlFailureTrack2');

tResult = Store.t;
xResult = Store.x;
strFileNameArr = sprintf('%s%s_timeAndPopArr.mat', folderName, fileName);
save(strFileNameArr, "tResult", 'xResult', 'eventBreakTimes')

end
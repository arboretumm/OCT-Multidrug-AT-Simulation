function [controlFailureTrack1, controlFailureTrack2] = sequentialDoubleBoundaryRunA_subPopSwitch(folderName, fileName, pStruct, populationNum, populationLim)
%file and folderName save the [ANYTHING ELSE?] and the figures generated in
%the filepath created 
close all;

%GENERAL RULES
    %use camelcase, not underscores for parameter names EXCEPTING things
    %that would have otherwise been sub or superscripted or need to split
    %apart chunks of a name
    %dear god comment you blessed idiot

    runType = 'sequentialDoubleBoundaryControlShapeA_subPopSwitch';
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
        [t, x] = ode45(@(t, x) controlACalcThreshold_controlBZero(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        eventBreakTimes(1) = Store.t(end);
        switchPoint1 = length(Store.t);
        %sprintf("phase 1 done")
        if t(end) < tSPAN(end)
            %phase 2 - control B at threshold and control A at max until
            %time-based switch
            %sprintf("phase 2 start")
            options2 = odeset('Events', @(t,x) eventSubpopulationLimit(t, x, p, populationNum, populationLim),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
            [t, x] = ode45(@(t, x) controlAMax_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options2);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
            %sprintf("phase 2 done")
            eventBreakTimes(2) = Store.t(end);
            switchPoint2 = length(Store.t);
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

%graphing concentration too
    concA = zeros(1, length(Store.t));
    concB = zeros(1, length(Store.t));

    concCalcA = ((p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4)) .* ((1-(p.kappa_threshold./p.kappa))))./(p.alphaA .* (Store.x(:,1) + Store.x(:,2)));
    concCalcB = (((p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4)) .* (1-(p.kappa_threshold./p.kappa))) - p.alphaA.*p.uMaxA.*(Store.x(:,1)+Store.x(:,2)))./(p.alphaB .* (Store.x(:,1) + Store.x(:,3)));


    for i=1:switchPoint1
        refA = concCalcA(i);
        if refA < 0
            concA(i) = 0;
        elseif refA > p.uMaxA
            concA(i) = p.uMaxB;
        else
            concA(i) = refA;
        end
    end
    
    
    
    for i=switchPoint1+1:switchPoint2
        refB = concCalcB(i);
        if refB < 0
            concB(i) = 0;
        elseif refB > p.uMaxB
            concB(i) = p.uMaxB;
        else
            concB(i) = refB;
        end
        concA(i) = p.uMaxA;
    end

    for i=switchPoint2+1:length(Store.t)
        concA(i) = p.uMaxA;
        concB(i) = p.uMaxB;
    end
    
fileName2 = sprintf('%s%s_concentrationGraph.fig', folderName, fileName);

gcf2 = figure('Visible', 'off');
hold on
    plot(Store.t, concA, 'LineWidth', 3)
    plot(Store.t, concB, 'LineWidth', 3)
    set(gcf2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')')
    set(gcf2, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
hold off

savefig(gcf2, fileName2);

% %parameter file save
% strFileName3 = sprintf('%s%s.mat', folderName, fileName);
% save(strFileName3, 'p', 'c', "runType");
[rowsz, colsz] = size(Store.x);
controlFailureTrack1 = 0;
controlFailureTrack2 = 0;
errorBreak = 1.001*(p.kappa_threshold);

cf1Index = find((Store.x(:, 4) > errorBreak), 1, 'first');
cf2Index = find((Store.x(:, 5) > errorBreak), 1, 'first');

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
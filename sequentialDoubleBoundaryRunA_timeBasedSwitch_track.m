function [controlFailureTrack1, controlFailureTrack2] = sequentialDoubleBoundaryRunA_timeBasedSwitch_track(folderName, fileName, pStruct, time1, time2)
%file and folderName save the [ANYTHING ELSE?] and the figures generated in
%the filepath created
tic
close all;

%GENERAL RULES
    %use camelcase, not underscores for parameter names EXCEPTING things
    %that would have otherwise been sub or superscripted or need to split
    %apart chunks of a name
    %dear god comment you blessed idiot

    runType = 'sequentialDoubleBoundaryControlShapeA_timeSwitch';
    %sequential double bound with A as the first boundary, B as 0 -> second
    %boundary, with switches between different phases time based

    %parameter definition
    %base equation build + instantiation of structure
        %structures now built outside of main run file (this one) and fed
        %in as full structures to this file in order to set up multi sim
        %runs easier

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
        options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        %sprintf("phase 1 start")
        [t, x] = ode45(@(t, x) controlACalcThreshold_controlBZero(t, x, p), [t(end) time1], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        %sprintf("phase 1 done")
        if t(end) < tSPAN(end)
            %phase 2 - control B at threshold and control A at max until
            %time-based switch
            %sprintf("phase 2 start")
            [t, x] = ode45(@(t, x) controlAMax_controlBCalcThreshold(t, x, p), [t(end) time2], x(end,:), options);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
            %sprintf("phase 2 done")
            if t(end) < tSPAN(end)
                %phase 3 - control B and control A at max until end of
                %simulation
                %sprintf("phase 3 start")
                [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
                %sprintf("phase 3 done")
            end
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

% gcf2 = figure(2);
% hold on;
% plot(Store.t, pStruct.concATrack, 'LineWidth', 3);
% plot(Store.t, pStruct.concBTrack, 'LineWidth', 3);
% axis([0 tSPAN(end) 0 inf]);
% title('Plot of drug concentrations', 'FontSize', 24);
% xlabel('time', 'FontSize', 18);
% legend('conc A', 'conc B');
% set(gcf2, 'units', 'points', 'position', [10, 10, 576, 576]); %points = 1/72", ~8" here?
% set(gca, 'FontSize', 16)
% hold off;
% 
% strFileName4 = sprintf('%s%s_conc.fig', folderName, fileName);
% strFileName5 = sprintf('%s%s_conc.png', folderName, fileName);
% saveas(figure(2),strFileName4); %saves population graph as fig
% saveas(figure(2),strFileName5); %saves population graph as png

% %parameter file save
% strFileName3 = sprintf('%s%s.mat', folderName, fileName);
% save(strFileName3, 'p', 'c', "runType");
% %parameter file save
% strFileName3 = sprintf('%s%s.mat', folderName, fileName);
% save(strFileName3, 'p', 'c', "runType");
[rowsz, colsz] = size(Store.x);
foundBreachC = false;
foundBreachD = false;
controlFailureTrack1 = 0;
controlFailureTrack2 = 0;

for q=1:rowsz
    if Store.x(q, 4) > p.kappa_threshold && foundBreachC == false
        disp('simulation broke threshold (4) at');
        disp(Store.t(q));
        controlFailureTrack1 = Store.t(q);
        foundBreachC = true;
    end
end

for q=1:rowsz
    if Store.x(q, 5) > p.kappa_threshold+1 && foundBreachD == false
        disp('simulation broke threshold (5) at');
        disp(Store.t(q));
        controlFailureTrack2 = Store.t(q);
        foundBreachD = true;
    end
end

%parameter file save with post calc
strFileName3 = sprintf('%s%s.mat', folderName, fileName);
save(strFileName3, 'p', 'runType', 'controlFailureTrack1');

toc

end
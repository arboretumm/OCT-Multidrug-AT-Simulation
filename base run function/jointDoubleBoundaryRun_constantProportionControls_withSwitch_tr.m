function [controlFailureTrack1, controlFailureTrack2] = jointDoubleBoundaryRun_constantProportionControls_withSwitch_try2(folderName, fileName, pStruct, constProp)
%file and folderName save the [ANYTHING ELSE?] and the figures generated in
%the filepath created
%tic
close all;

%GENERAL RULES
    %use camelcase, not underscores for parameter names EXCEPTING things
    %that would have otherwise been sub or superscripted or need to split
    %apart chunks of a name
    %dear god comment you blessed idiot

    runType = 'jointDoubleBoundaryRun_constantProportionControls';
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
    %t=tSPAN(1);
    x=initCond;

    options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
    [t, x] = ode45(@(t, x) jointBoundaryEqn_constantProportionControls_withSwitch_try2(t, x, p, constProp), tSPAN, x(end,:), options);
    Store.t=[Store.t ;t(2:end)];
    Store.x=[Store.x; x(2:end,:)];


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
    saveas(gcf1,strFileName1); %saves population graph as fig
    saveas(gcf1,strFileName2); %saves population graph as png

    %build concentration vectors
    


    %CONCENTRATION GRAPHING
    gcf2 = figure('Visible', 'off');
    concentrationA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4)))./(p.alphaA.*(Store.x(:,1)+Store.x(:,2))+p.alphaB.*constProp.*(Store.x(:,1)+Store.x(:,3)));
    concentrationB = constProp*concentrationA;
    hold on;
    plot(Store.t, concentrationA, 'LineWidth', 2.5);
    plot(Store.t, concentrationB, 'LineWidth', 2.5);
    title('Plot of drug concentration', 'FontSize', 24);
    xlabel('time (hours)', 'FontSize', 18);
    ylabel('drug concentration (nM)', 'FontSize', 18);
    legend('Drug A', 'Drug B');
    set(gcf2, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    axis([0 tSPAN(end) 0 100]);
    set(gcf2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    hold off;

    strFileName5 = sprintf('%s%s%s.fig', folderName, fileName, "_drugConc");
    strFileName6 = sprintf('%s%s%s.png', folderName, fileName, "_drugConc");
    saveas(gcf2,strFileName5); %saves population graph as fig
    saveas(gcf2,strFileName6); %saves population graph as png

    %CONCENTRATION GRAPHING
    gcf3 = figure('Visible', 'off');
    hold on;
    plot(Store.t, concentrationA, 'LineWidth', 2.5);
    plot(Store.t, concentrationB, 'LineWidth', 2.5);
    title('Plot of drug concentration, zoomed', 'FontSize', 24);
    xlabel('time (hours)', 'FontSize', 18);
    ylabel('drug concentration (nM)', 'FontSize', 18);
    legend('Drug A', 'Drug B');
    set(gcf3, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    axis([0 tSPAN(end) 0 1]);
    set(gcf3, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    hold off;

    strFileName7 = sprintf('%s%s%s.fig', folderName, fileName, "_drugConcZoom");
    strFileName8 = sprintf('%s%s%s.png', folderName, fileName, "_drugConcZoom");
    saveas(gcf3,strFileName7); %saves population graph as fig
    saveas(gcf3,strFileName8); %saves population graph as png
%find breach by pop 4
[rowsz, colsz] = size(Store.x);
foundBreachC = false;
foundBreachD = false;
controlFailureTrack1 = 0;
controlFailureTrack2 = 0;

for q=1:rowsz
    if Store.x(q, 4) > p.kappa_threshold && foundBreachC == false
        %disp('simulation broke threshold (4) at');
        %disp(Store.t(q));
        controlFailureTrack1 = Store.t(q);
        foundBreachC = true;
    end
end
for q=1:rowsz
    if Store.x(q, 5) > p.kappa_threshold+1 && foundBreachD == false
        %disp('simulation broke threshold (5) at');
        %disp(Store.t(q));
        controlFailureTrack2 = Store.t(q);
        foundBreachD = true;
    end
end

%add graph of concentration over time? 

%parameter file save with post calc
strFileName3 = sprintf('%s%s.mat', folderName, fileName);
save(strFileName3, 'p', 'runType', 'controlFailureTrack1');

%toc
end
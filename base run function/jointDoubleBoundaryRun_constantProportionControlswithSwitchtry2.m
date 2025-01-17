function [controlFailureTrack1, controlFailureTrack2] = jointDoubleBoundaryRun_constantProportionControlswithSwitchtry2(folderName, fileName, pStruct, constProp)
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
    concentrationAStart = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4)))./(p.alphaA.*(Store.x(:,1)+Store.x(:,2))+p.alphaB.*constProp.*(Store.x(:,1)+Store.x(:,3)));
    concentrationBStart = constProp.*concentrationAStart;

    concentrationA = zeros(1, length(concentrationAStart));
    concentrationB = zeros(1, length(concentrationBStart));

    for i=1:length(concentrationAStart)
        concA = concentrationAStart(i);
        concB = concentrationBStart(i);
        if concA < 0
            concentrationA(1, i) = 0;
            concA = 0;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(i,1) + p.lambda_m .* Store.x(i,2) + p.lambda_s .* Store.x(i,3) + p.lambda_d .* Store.x(i,4)) - (p.alphaA.*concA.*(Store.x(i,1)+Store.x(i,2))))./(p.alphaB.*(Store.x(i, 1)+Store.x(i, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concentrationB(1, i) = 0;
            elseif recalcB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
            else
                concentrationB(1, i) = recalcB;
            end
        elseif concA > p.uMaxA
            concentrationA(1, i) = p.uMaxA;
            concA = p.uMaxA;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(i,1) + p.lambda_m .* Store.x(i,2) + p.lambda_s .* Store.x(i,3) + p.lambda_d .* Store.x(i,4)) - (p.alphaA.*concA.*(Store.x(i,1)+Store.x(i,2))))./(p.alphaB.*(Store.x(i, 1)+Store.x(i, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concentrationB(1, i) = 0;
            elseif recalcB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
            else
                concentrationB(1, i) = recalcB;
            end
        else
            if concB < 0
                concentrationB(1, i) = 0;
                concB = 0;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(i,1) + p.lambda_m .* Store.x(i,2) + p.lambda_s .* Store.x(i,3) + p.lambda_d .* Store.x(i,4)) - (p.alphaB.*concB.*(Store.x(i,1)+Store.x(i,3))))./(p.alphaA.*(Store.x(i,1)+Store.x(i,2)));
                if recalcA < 0
                    concentrationA(1, i) = 0;
                elseif recalcA > p.uMaxA
                    concentrationA(1, i) = p.uMaxA;
                else
                    concentrationA(1, i) = recalcA;
                end
            elseif concB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
                concB = p.uMaxB;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(i,1) + p.lambda_m .* Store.x(i,2) + p.lambda_s .* Store.x(i,3) + p.lambda_d .* Store.x(i,4)) - (p.alphaB.*concB.*(Store.x(i,1)+Store.x(i,3))))./(p.alphaA.*(Store.x(i,1)+Store.x(i,2)));
                if recalcA < 0
                    concentrationA(1, i) = 0;
                elseif recalcA > p.uMaxA
                    concentrationA(1, i) = p.uMaxA;
                else
                    concentrationA(1, i) = recalcA;
                end
            else
                concentrationA(1, i) = concA;
                concentrationB(1, i) = concB;
            end
        end
    end

    %CONCENTRATION GRAPHING
    gcf2 = figure('Visible', 'off');
%     concentrationA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4)))./(p.alphaA.*(Store.x(:,1)+Store.x(:,2))+p.alphaB.*constProp.*(Store.x(:,1)+Store.x(:,3)));
%     concentrationB = constProp*concentrationA;
    hold on;
    plot(Store.t, concentrationA, 'LineWidth', 2.5);
    plot(Store.t, concentrationB, 'LineWidth', 2.5);
    title('Plot of drug concentration', 'FontSize', 24);
    xlabel('time (hours)', 'FontSize', 18);
    ylabel('drug concentration (nM)', 'FontSize', 18);
    legend('Drug A', 'Drug B');
    set(gcf2, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    axis([0 tSPAN(end) 0 inf]);
    set(gcf2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    hold off;

    strFileName5 = sprintf('%s%s%s.fig', folderName, fileName, "_drugConc");
    strFileName6 = sprintf('%s%s%s.png', folderName, fileName, "_drugConc");
    saveas(gcf2,strFileName5); %saves population graph as fig
    saveas(gcf2,strFileName6); %saves population graph as png

    % %CONCENTRATION GRAPHING ON TOP OF THE OTHER GRAPHING
    % gcf4 = figure('Visible', 'off');
    % hold on;
    % colororder({'k','k'})
    % yyaxis left;
    % plot(Store.t, Store.x(:,1), 'LineWidth', 2, 'Color',"#0072BD", 'LineStyle',"-", "Marker", 'none');
    % plot(Store.t, Store.x(:,2), 'LineWidth', 2, 'Color',"#D95319", 'LineStyle',"-", "Marker", 'none');
    % plot(Store.t, Store.x(:,3), 'LineWidth', 2, 'Color',"#EDB120", 'LineStyle',"-", "Marker", 'none'); 
    % plot(Store.t, Store.x(:,4), 'LineWidth', 2, 'Color',"#7E2F8E", 'LineStyle',"-", "Marker", 'none');
    % plot(Store.t, Store.x(:,5), 'LineWidth', 2, 'Color',"#77AC30", 'LineStyle',"-", "Marker", 'none');
    % axis([0 tSPAN(end) 0 inf]);
    % title('Cell populations and drug concentration', 'FontSize', 24);
    % xlabel('time (hours)', 'FontSize', 18);
    % ylabel('cell number', 'FontSize', 18);
    % 
    % yyaxis right;
    % plot(Store.t, concentrationA, 'LineWidth', 2.5, 'LineStyle','--');
    % plot(Store.t, concentrationB, 'LineWidth', 2.5, 'LineStyle',':');
    % ylabel('drug concentration (nM)', 'FontSize', 18);
    % 
    % set(gcf4, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    % set(gca, 'FontSize', 16)
    % axis([0 tSPAN(end) 0 1]);
    % set(gcf4, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    % hold off;
    % strFileName9 = sprintf('%s%s%s.fig', folderName, fileName, "_drugConcAndPopulation");
    % strFileName10 = sprintf('%s%s%s.png', folderName, fileName, "_drugConcAndPopulation");
    % saveas(gcf4,strFileName9); %saves population graph as fig
    % saveas(gcf4,strFileName10); %saves population graph as png
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

tResult = Store.t;
xResult = Store.x;
strFileNameArr = sprintf('%s%s_timeAndPopArr.mat', folderName, fileName);
save(strFileNameArr, "tResult", 'xResult')
%toc
end
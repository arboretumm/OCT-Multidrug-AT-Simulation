function [controlFailureTrack1, controlFailureTrack2] = sequentialDoubleBoundaryRunB_concDrugSwitch_withAOption(folderName, fileName, pStruct, drugConc1, drugConc2)
%file and folderName save the [ANYTHING ELSE?] and the figures generated in
%the filepath created
%tic
close all;

%GENERAL RULES
    %use camelcase, not underscores for parameter names EXCEPTING things
    %that would have otherwise been sub or superscripted or need to split
    %apart chunks of a name
    %dear god comment you blessed idiot

    runType = 'sequentialDoubleBoundaryControlShapeA_concDrugSwitch';
    %sequential double bound with A as the first boundary, B as 0 -> second
    %boundary, with switches between different phases based on the
    %concentration of the boundary drug primarily

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
    concentrationA = [];
    concentrationB = [];
    initialDrugAVal = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC)- 0.*p.alphaB.*(p.naivePopIC + p.mutBPopIC))./(p.alphaA .* (p.naivePopIC + p.mutAPopIC));
    initialDrugBVal = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC)- 0.*p.alphaA.*(p.naivePopIC + p.mutAPopIC))./(p.alphaB .* (p.naivePopIC + p.mutBPopIC));
    
    if initialDrugBVal < drugConc2
        if t(end) < tSPAN(end)
        %phase 1 - control B at 0 and control A at threshold until
        %time-based switch
        options = odeset('Events', @(t,x) eventGoalB_ASub(t, x, p, drugConc2, 0),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
        %sprintf("phase 1 start")
        [t, x] = ode45(@(t, x) controlAZero_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        lengthVal = length(Store.t);
        concentrationAHold = zeros(1, lengthVal);
        concentrationA = [concentrationA; concentrationAHold];
        concentrationBHold = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(:, 1) + p.lambda_m .* Store.x(:, 2) + p.lambda_s .* Store.x(:, 3) + p.lambda_d .* Store.x(:, 4))- 0.*p.alphaA.*(Store.x(:, 1) + Store.x(:, 2)))./(p.alphaB .* (Store.x(:, 1) + Store.x(:, 3)));
        concentrationB = [concentrationB; concentrationBHold'];
        %sprintf("phase 1 done")
            if t(end) < tSPAN(end)
                %phase 2 - control B at threshold and control A at max until
                %time-based switch
                %sprintf("phase 2 start")
                options2 = odeset('Events', @(t, x) eventGoalA_BSub(t, x, p, drugConc1, p.uMaxB), 'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                [t, x] = ode45(@(t, x) controlACalcThreshold_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options2);
                xHold2 = x(2:end,:);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
                lengthVal = length(t) - 1;
                concentrationAHold=((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xHold2(:,1) + p.lambda_m .* xHold2(:,2) + p.lambda_s .* xHold2(:,3) + p.lambda_d .* xHold2(:,4))- p.uMaxB.*p.alphaB.*(xHold2(:,1) + xHold2(:,3)))./(p.alphaA .* (xHold2(:,1) + xHold2(:,2)));
                concentrationA = [concentrationA, concentrationAHold'];
                concentrationBHold2 = p.uMaxB*ones(1, lengthVal);
                concentrationB = [concentrationB, concentrationBHold2];
                %sprintf("phase 2 done")
                if t(end) < tSPAN(end)
                    %phase 3 - control B and control A at max until end of
                    %simulation
                    %sprintf("phase 3 start")
                    options3 = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                    [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options3);
                    Store.t=[Store.t ;t(2:end)];
                    Store.x=[Store.x; x(2:end,:)];
                    lengthVal3 = length(t) - 1;
                    concentrationAHold3=p.uMaxA.*ones(1,lengthVal3);
                    concentrationA = [concentrationA, concentrationAHold3];
                    concentrationBHold3=p.uMaxB*ones(1,lengthVal3);
                    concentrationB = [concentrationB, concentrationBHold3];
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
        saveas(gcf1,strFileName1); %saves population graph as fig
        saveas(gcf1,strFileName2); %saves population graph as png
        
        %CONCENTRATION GRAPHING
        gcf2 = figure('Visible', 'off');
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
        %         disp('simulation broke threshold (4) at');
        %         disp(Store.t(q));
                controlFailureTrack1 = Store.t(q);
                foundBreachC = true;
            end
        end
        
        for q=1:rowsz
            if Store.x(q, 5) > p.kappa_threshold+1 && foundBreachD == false
        %         disp('simulation broke threshold (5) at');
        %         disp(Store.t(q));
                controlFailureTrack2 = Store.t(q);
                foundBreachD = true;
            end
        end
        
        %parameter file save with post calc
        strFileName3 = sprintf('%s%s.mat', folderName, fileName);
        save(strFileName3, 'p', 'runType', 'controlFailureTrack1');

    elseif initialDrugAVal < drugConc1 && initialDrugBVal >= drugConc2
        if t(end) < tSPAN(end)
        %phase 1 - control B at 0 and control A at threshold until
        %time-based switch
        options = odeset('Events', @(t,x) eventGoalA_BSub(t, x, p, drugConc1, 0),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
        %sprintf("phase 1 start")
        [t, x] = ode45(@(t, x) controlACalcThreshold_controlBZero(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        lengthVal = length(Store.t);
        concentrationAHold = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* Store.x(:,1) + p.lambda_m .* Store.x(:,2) + p.lambda_s .* Store.x(:,3) + p.lambda_d .* Store.x(:,4))- 0.*p.alphaB.*(Store.x(:,1) + Store.x(:,3)))./(p.alphaA .* (Store.x(:,1) + Store.x(:,2)));
        concentrationA = [concentrationA; concentrationAHold'];
        concentrationBHold = zeros(1, lengthVal);
        concentrationB = [concentrationB; concentrationBHold];
            if t(end) < tSPAN(end)
                %phase 2 - control B at threshold and control A at max until
                %time-based switch
                %sprintf("phase 2 start")
                options2 = odeset('Events', @(t, x) eventGoalB_ASub(t, x, p, drugConc1, p.uMaxB), 'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                [t, x] = ode45(@(t, x) controlAMax_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options2);
                xHold2 = x(2:end,:);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
                lengthVal = length(t) - 1;
                concentrationAHold=p.uMaxA*ones(1, lengthVal);
                concentrationA = [concentrationA, concentrationAHold];
                %CONCENTRATIONBHOLD IS GENERATING A SQUARE MATRIX FOR SOME
                %REASON
                concentrationBHold=((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .*xHold2(:, 1) + p.lambda_m .* xHold2(:, 2) + p.lambda_s .* xHold2(:, 3) + p.lambda_d .* xHold2(:, 4))- p.uMaxA.*p.alphaA.*(xHold2(:, 1) + xHold2(:, 2)))./(p.alphaB .* (xHold2(:, 1) + xHold2(:, 3)));
                concentrationB = [concentrationB, concentrationBHold'];
                %sprintf("phase 2 done")
                if t(end) < tSPAN(end)
                    %phase 3 - control B and control A at max until end of
                    %simulation
                    %sprintf("phase 3 start")
                    options3 = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                    [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options3);
                    Store.t=[Store.t ;t(2:end)];
                    Store.x=[Store.x; x(2:end,:)];
                    lengthVal3 = length(t) - 1;
                    concentrationAHold3=p.uMaxA.*ones(1,lengthVal3);
                    concentrationA = [concentrationA, concentrationAHold3];
                    concentrationBHold3=p.uMaxB*ones(1,lengthVal3);
                    concentrationB = [concentrationB, concentrationBHold3];
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
        saveas(gcf1,strFileName1); %saves population graph as fig
        saveas(gcf1,strFileName2); %saves population graph as png
        
        %CONCENTRATION GRAPHING
        gcf2 = figure('Visible', 'off');
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
        %         disp('simulation broke threshold (4) at');
        %         disp(Store.t(q));
                controlFailureTrack1 = Store.t(q);
                foundBreachC = true;
            end
        end
        
        for q=1:rowsz
            if Store.x(q, 5) > p.kappa_threshold+1 && foundBreachD == false
        %         disp('simulation broke threshold (5) at');
        %         disp(Store.t(q));
                controlFailureTrack2 = Store.t(q);
                foundBreachD = true;
            end
        end
        
        %parameter file save with post calc
        strFileName3 = sprintf('%s%s.mat', folderName, fileName);
        save(strFileName3, 'p', 'runType', 'controlFailureTrack1');

    else
        controlFailureTrack1 = -1000;
        controlFailureTrack2 = -1000;
    end
end
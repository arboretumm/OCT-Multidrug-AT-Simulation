function [thresholdBreak] = maximalTreatmentSequential(popStructFile, runFile, saveFile, drugConcA, drugConcB, timestep)
        load(popStructFile, "popStruct"); %population structure, required for equation functions usually
        %establishes the time and population of the actual original run,
        %will be used to back calculate the concentration of the original
        %run
        load(runFile, "tResult");
        load(runFile, "xResult");

        p = popStruct;
    
        endTime = p.cellTime + 1;

        timeVector = 1:timestep:endTime;
        concAVec = drugConcA.*ones(length(timeVector));
        concBVec = drugConcB.*ones(length(timeVector));
        concZeroVec = zeros(length(timeVector));

        %feed new concentration vectors into the equations, then compare
        %previous simulation with new simulation & fail time

        initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];
    
        %solving/graphing system given parameters
        tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
        %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
    
        Store.t=[];
        Store.x=[];
        %t=tSPAN(1);
        x=initCond;
        x2 = initCond; 
        xMin = initCond(5);

    %find minimum value
        optionsFind = odeset('Events', @(t, x) eventFindTreatmentMinima_PracticalityA(t, x, p, timeVector, concAVec, concZeroVec), 'RelTol', 1e-4, 'NonNegative', [1 2 3 4 5], 'Refine', 10);
        [t2, x2, te_a, xe_a, ie_a] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concAVec, concZeroVec), tSPAN, x(end,:), optionsFind);
        if length(te_a) ~= 0
            if xMin > xe_a(5)
                xMin = xe_a(5);
            end
        end
    %maintain regimen until 20% increase for RECIST
        options = odeset('Events', @(t,x) eventSwitchMaxRECIST(t, x, p, xMin),'RelTol', 1e-4, 'NonNegative', [1 2 3 4 5], 'Refine', 10);
        [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concAVec, concZeroVec), tSPAN, x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];
        if t(end) < tSPAN(end)
            options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
            [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concZeroVec, concBVec), [t(end) tSPAN(end)], x(end,:), options);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
        end

        %GRAPHING

            gcf2 = figure('Visible', 'off');
            hold on;
            line_color = ['#0072BD', 	"#D95319", 	"#EDB120", 	"#7E2F8E", 	"#77AC30", "#9fc8e2", 	"#d8b7a3", 	"#c4a55c", 	"#c9b4cd", 	"#bfcdac"];
            plot(Store.t, Store.x, 'LineWidth', 2.5);
            plot(tResult, xResult, 'LineWidth', 2.5, 'LineStyle', '--')
            colororder(line_color)
            axis([0 tSPAN(end) 0 inf]);
            title('Plot of cell populations - Adjusted Regimen', 'FontSize', 24);
            xlabel('time (hours)', 'FontSize', 18);
            ylabel('cell number', 'FontSize', 18);
            legend('susceptible', 'B resist', 'A resist', 'both resist', 'total');
            set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
            set(gca, 'FontSize', 16)
            set(gcf, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
            hold off
            %save figure
            strFileName3 = sprintf('%s%s.fig', saveFile, "_cellPopulationApproxTx");
            strFileName4 = sprintf('%s%s.png', saveFile, "_cellPopulationApproxTx");
            saveas(gcf2,strFileName3); %saves population graph as fig
            saveas(gcf2,strFileName4); %saves population graph as png

       %find breach by pop 4
            [rowsz, colsz] = size(Store.x);
            foundBreachC = false;
            foundBreachD = false;
            controlFailureTrack1 = 0;
            controlFailureTrack2 = 0;
            errorBreak = 1.001*(p.kappa_threshold);
            
            cf1Index = find((Store.x(:, 4) > errorBreak), 1, 'first');
            cf2Index = find((Store.x(:, 5) > errorBreak), 1, 'first');

            controlFailureTrack1 = Store.t(cf1Index);
            controlFailureTrack2 = Store.t(cf2Index);

            thresholdBreak = controlFailureTrack2;
end 
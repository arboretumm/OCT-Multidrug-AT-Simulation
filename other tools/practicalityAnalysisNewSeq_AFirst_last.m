function [thresholdBreak] = practicalityAnalysisNewSeq_AFirst_last(popStructFile, runFile, saveFile, timestep, dosingIncrement)
    %establish run we're comparing to
        load(popStructFile, "popStruct"); %population structure, required for equation functions usually
        %establishes the time and population of the actual original run,
        %will be used to back calculate the concentration of the original
        %run
        p = popStruct;
        load(runFile, "tResult");
        load(runFile, "xResult");
        load(runFile, "eventBreakTimes");

    %back calculate the concentration of the original run
        
       [rows, columns] = size(xResult);

        concentrationA = zeros(1, rows);
        concentrationB = zeros(1, rows);
        
        breakOne = eventBreakTimes(1);
        breakTwo = eventBreakTimes(2);

        timeVal = tResult(1);
        timeIndex = 1;
        initialDrugAVal = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC)- 0.*p.alphaB.*(p.naivePopIC + p.mutBPopIC))./(p.alphaA .* (p.naivePopIC + p.mutAPopIC));
        initialDrugBVal = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC)- 0.*p.alphaA.*(p.naivePopIC + p.mutAPopIC))./(p.alphaB .* (p.naivePopIC + p.mutBPopIC));

        while timeVal <= breakOne
            concentrationA(timeIndex) = ((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* ((1-(p.kappa_threshold./p.kappa))))/(p.alphaA .* (xResult(timeIndex, 1) + xResult(timeIndex, 2)));
            concentrationB(timeIndex) = 0;

            timeIndex = timeIndex + 1;
            timeVal = tResult(timeIndex);
        end

        while timeVal <= breakTwo
            concentrationB(timeIndex) = (((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* (1-(p.kappa_threshold./p.kappa))) - (p.alphaA.*p.uMaxA.*(xResult(timeIndex, 1)+xResult(timeIndex, 2))))/(p.alphaB .* (xResult(timeIndex, 1) + xResult(timeIndex, 3)));
            concentrationA(timeIndex) = p.uMaxA;

            timeIndex = timeIndex + 1;
            timeVal = tResult(timeIndex);
        end

        for i=timeIndex:length(tResult)
            concentrationA(i) = p.uMaxA;
            concentrationB(i) = p.uMaxB;
        end
        % 
        % 
        % while timeVal <= breakOne
        %     if initialDrugBVal < p.uMaxB
        %         concentrationB(timeIndex) = ((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* ((1-(p.kappa_threshold./p.kappa))))/(p.alphaB .* (xResult(timeIndex, 1) + xResult(timeIndex, 3)));
        %         concentrationA(timeIndex) = 0;
        %     else
        %         concentrationA(timeIndex) = ((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* ((1-(p.kappa_threshold./p.kappa))))/(p.alphaA .* (xResult(timeIndex, 1) + xResult(timeIndex, 2)));
        %         concentrationB(timeIndex) = 0;
        %     end
        % 
        %     timeIndex = timeIndex + 1;
        %     timeVal = tResult(timeIndex);
        % end
        % 
        % while timeVal <= breakTwo
        %     if initialDrugAVal < p.uMaxB
        %         concentrationA(timeIndex) = (((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* ((1-(p.kappa_threshold./p.kappa)))) - (p.alphaB.*p.uMaxB.*(xResult(timeIndex, 1) + xResult(timeIndex, 3))))/(p.alphaA .* (xResult(timeIndex, 1) + xResult(timeIndex, 2)));
        %         concentrationB(timeIndex) = p.uMaxB;
        %     else
        %         concentrationB(timeIndex) = (((p.lambda_n .* xResult(timeIndex, 1) + p.lambda_m .* xResult(timeIndex, 2) + p.lambda_s .* xResult(timeIndex, 3) + p.lambda_d .* xResult(timeIndex, 4)) .* (1-(p.kappa_threshold./p.kappa))) - (p.alphaA.*p.uMaxA.*(xResult(timeIndex, 1)+xResult(timeIndex, 2))))/(p.alphaB .* (xResult(timeIndex, 1) + xResult(timeIndex, 3)));
        %         concentrationA(timeIndex) = p.uMaxA;
        %     end
        % 
        %     timeIndex = timeIndex + 1;
        %     timeVal = tResult(timeIndex);
        % end
        % 
        % for i=timeIndex:length(tResult)
        %     concentrationA(i) = p.uMaxA;
        %     concentrationB(i) = p.uMaxB;
        % end
        
    %use timestep & interpolation to identify concentration for our
    %stepwise treatment regimen -- play w/floor v ceiling rounding?
        timeChunks = floor(p.cellTime./timestep);
        endTime = timeChunks*timestep;

        timeVals = 0:timestep:endTime;

        interpConcA = interp1(tResult, concentrationA, timeVals, 'makima');
        interpConcB = interp1(tResult, concentrationB, timeVals, 'makima');

        acc = dosingIncrement;

        roundedConcA = ceil(interpConcA/acc)*acc;
        roundedConcB = ceil(interpConcB/acc)*acc;
        
        for k =1:length(roundedConcA)
            if roundedConcA(k) > p.uMaxA
                roundedConcA(k) = p.uMaxA;
            elseif roundedConcA(k) < 0
                roundedConcA(k) = 0;
            end
        end

        for m =1:length(roundedConcB)
            if roundedConcB(m) > p.uMaxB
                roundedConcB(m) = p.uMaxB;
            elseif roundedConcB(k) < 0
                roundedConcB(k) = 0;
            end
        end
        % 

            gcf1 = figure('Visible', 'off');
            hold on
            plot(tResult, concentrationA, 'LineWidth', 3, 'LineStyle', '-', 'Color', '#7441d3')
            plot(timeVals, roundedConcA, 'LineWidth', 3, 'LineStyle', '-.', 'Color', '#00a1bd')

            plot(tResult, concentrationB, 'LineWidth', 3, 'LineStyle', '-', 'Color', '#77AC30')
            plot(timeVals, roundedConcB, 'LineWidth', 3, 'LineStyle', '-.', 'Color', '#d6b606')
            title('Plot of drug concentrations', 'FontSize', 24);
            xlabel('time (hours)', 'FontSize', 18);
            ylabel('concentration (nM)', 'FontSize', 18);
            legend('[drug A]', 'practical [drug A]', '[drug B]', 'practical [drug B]');
            set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
            set(gca, 'FontSize', 16)
            set(gcf, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
            hold off
            %save figure
            strFileName1 = sprintf('%s%s.fig', saveFile, "_drugConcentrationApprox");
            strFileName2 = sprintf('%s%s.png', saveFile, "_drugConcentrationApprox");
            saveas(gcf1,strFileName1); %saves population graph as fig
            saveas(gcf1,strFileName2); %saves population graph as png

            disp("defined all our sequential variables")
     
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
    
        options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVals, roundedConcA, roundedConcB), tSPAN, x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];

        disp("made it through the ODE")

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
            
            % cf1Index = find((Store.x(:, 4) > errorBreak), 1, 'first');
            % cf2Index = find((Store.x(:, 5) > errorBreak), 1, 'first');
            
            cf1Index = find((Store.x(:, 4) <= errorBreak), 1, 'last');
            cf2Index = find((Store.x(:, 5) <= errorBreak), 1, 'last');

            controlFailureTrack1 = Store.t(cf1Index);
            controlFailureTrack2 = Store.t(cf2Index);
            
            thresholdBreak = controlFailureTrack2;
        
end
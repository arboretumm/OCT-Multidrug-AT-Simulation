function [thresholdBreak] = practicalityAnalysisJoint(popStructFile, runFile, saveFile, timestep, constProp, dosingIncrement)
    %establish run we're comparing to
        load(popStructFile, "popStruct"); %population structure, required for equation functions usually
        %establishes the time and population of the actual original run,
        %will be used to back calculate the concentration of the original
        %run
        p = popStruct;
        load(runFile, "tResult");
        load(runFile, "xResult")

    %back calculate the concentration of the original run -- 
        concentrationAStart = ((1 - p.kappa_threshold/p.kappa)*(p.lambda_n .* xResult(:,1) + p.lambda_m .* xResult(:,2) + p.lambda_s .* xResult(:,3) + p.lambda_d .* xResult(:,4)))./(p.alphaA.*(xResult(:,1)+xResult(:,2))+p.alphaB.*constProp.*(xResult(:,1)+xResult(:,3)));
        concentrationBStart = constProp*concentrationAStart;

        concentrationA = zeros(1, length(concentrationAStart));
        concentrationB = zeros(1, length(concentrationBStart));

    for i=1:length(concentrationAStart)
        concA = concentrationAStart(i);
        concB = concentrationBStart(i);
        if concA < 0
            concentrationA(1, i) = 0;
            concA = 0;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaA.*concA.*(xResult(i,1)+xResult(i,2))))./(p.alphaB.*(xResult(i, 1)+xResult(i, 3))); %insertequationusingconcAinstead
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
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaA.*concA.*(xResult(i,1)+xResult(i,2))))./(p.alphaB.*(xResult(i, 1)+xResult(i, 3))); %insertequationusingconcAinstead
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
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaB.*concB.*(xResult(i,1)+xResult(i,3))))./(p.alphaA.*(xResult(i,1)+xResult(i,2)));
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
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaB.*concB.*(xResult(i,1)+xResult(i,3))))./(p.alphaA.*(xResult(i,1)+xResult(i,2)));
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
            end
        end

        for m =1:length(roundedConcB)
            if roundedConcB(m) > p.uMaxB
                roundedConcB(m) = p.uMaxB;
            end
        end
        
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
            set(gcf1, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
            hold off
            %save figure
            strFileName1 = sprintf('%s%s.fig', saveFile, "_drugConcentrationApprox");
            strFileName2 = sprintf('%s%s.png', saveFile, "_drugConcentrationApprox");
            saveas(gcf1,strFileName1); %saves population graph as fig
            saveas(gcf1,strFileName2); %saves population graph as png
     
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

            set(gcf2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
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
            
            for q=1:rowsz
                if Store.x(q, 4) > errorBreak && foundBreachC == false
                    %disp('simulation broke threshold (4) at');
                    %disp(Store.t(q));
                    controlFailureTrack1 = Store.t(q);
                    foundBreachC = true;
                end
            end
            for q=1:rowsz
                if Store.x(q, 5) > errorBreak && foundBreachD == false
                    %disp('simulation broke threshold (5) at');
                    %disp(Store.t(q));
                    controlFailureTrack2 = Store.t(q);
                    foundBreachD = true;
                end
            end

            thresholdBreak = controlFailureTrack2;
        
end
function [failureTime, switchpoint] = multirun_maximalOCTShape_Bfirst(uMax, kR, folderBase)
    %population params that are staying the same

        %equation parameters
                %drug effect -- rn osimertinib and afatinib
                alphaA = 0.00155; %0.06; nM/hr
                alphaB = 0.00155; %0.06; nM/hr
            %population growth rate
                gammaN = 0.031; %0.031; %only EGFR+
                gammaA = 0.031;
                gammaB = 0.031;
                gammaD = 0.031; %0.011; %EGFR/T790M/C797S
                %no current difference between EGFR+ from L858R and ex19del
            %carrying capacity + population parameters
                kappa = 40000; %actual carrying capacity of cell culture // 40k 96MW
                thresholdRatio = kR;
                kappa_threshold = thresholdRatio.*kappa; %3500000; %goal capacity
            %drug limits
                uMaxA = uMax; %nM
                uMaxB = uMax; %nM
            %time parameters
                cellTime = 2400; %number of timesteps eqns are running for (in hours?)
            %initial conditions
                naivePopIC = 0.685.*kappa_threshold; 
                jointMutPopIC = 0.3.*kappa_threshold; 
                doubleMutPopIC = 0.015.*kappa_threshold;

                folderName = sprintf("%sotherOCTshapes/umax_%s_thresh_%s/", folderBase, num2str(uMax), num2str(kR));

            
    %ratios for looping through
        ratio_2naught_3naught = [.05, .25, .5, .75, .95];
        %ratio_2naught_3naught = [.5];

    endTime = cellTime + 1;
    timestep = 10;
    timeVector = 1:timestep:endTime;

    bestThresholdTimeArr = [];
    bestSwitchTimepointArr = [];

    for q = 1:length(ratio_2naught_3naught)
        %establish population
            initCondA = (1 - ratio_2naught_3naught(q)).*jointMutPopIC;
            initCondB = ratio_2naught_3naught(q).*jointMutPopIC;

         %pop struct & file names
            %rename file structure check on
            % fileName = sprintf("individualGraphs/test_x2x3ratio%s_gammaEqual", num2str(ratio_2naught_3naught(q)));
            % populationFileName = sprintf("%sindividualGraphs/population_x2x3ratio%s_gammaEqual.mat", folderName, num2str(ratio_2naught_3naught(q)));
            p = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, gammaN, gammaA, gammaB, gammaD, kappa, kappa_threshold, cellTime, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
            % save(populationFileName, "popStruct");

            concBVec = uMaxB.*ones(1, length(timeVector));
            concAVec = zeros(1, length(timeVector));
            
        
            failureVec = zeros(1, length(timeVector));

        for i = 1:length(timeVector)
            concAVecAlt = concAVec;
            concAVecAlt(i:end) = uMaxA;
    
            %feed new concentration vectors into the equations, then compare
            %previous simulation with new simulation & fail time
    
            initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];
        
            %solving/graphing system given parameters
            tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
            %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        
            Store.t=[];
            Store.x=[];
            x=initCond;
        
            options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
            [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concAVecAlt, concBVec), tSPAN, x(end,:), options);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];
    
    
    %         %find breach by pop 4 (last threshold break version)
    %         [rowsz, colsz] = size(Store.x);
    %         foundBreachC = false;
    %         foundBreachD = false;
    %         controlFailureTrack1 = 0;
    %         controlFailureTrack2 = 0;
    %         errorBreak = 1.001*(p.kappa_threshold);
    %         
    %         cf1Index = find((Store.x(:, 4) <= errorBreak), 1, 'last');
    %         cf2Index = find((Store.x(:, 5) <= errorBreak), 1, 'last');
    % 
    %         controlFailureTrack1 = Store.t(cf1Index);
    %         controlFailureTrack2 = Store.t(cf2Index);
    % 
    %         thresholdBreak = controlFailureTrack2;
    
            %find breach by pop 4 (first threshold break version)
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
    
            failureVec(i) = thresholdBreak;
        end
        [bestThresholdBreak, ind] = max(failureVec);
        bestSwitchTime = timeVector(ind);

        figureName = sprintf('%sindividualMaxMaxB/test_x2x3ratio%s_gammaEqual.fig', folderName, num2str(ratio_2naught_3naught(q)));
        
        concAVecBest = concAVec;
        concAVecBest(ind:end) = uMaxA;

        initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];
        
        %solving/graphing system given parameters
        tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
        %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
    
        Store.t=[];
        Store.x=[];
        x=initCond;
    
        options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concAVecBest, concBVec), tSPAN, x(end,:), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];

        
            gcf1 = figure('Visible', 'off');
            plot(Store.t, Store.x, 'LineWidth', 2.5);
            axis([0 tSPAN(end) 0 inf]);
            title('Plot of cell populations - Max/Max A First', 'FontSize', 24);
            xlabel('time (hours)', 'FontSize', 18);
            ylabel('cell number', 'FontSize', 18);
            legend('susceptible', 'B resist', 'A resist', 'both resist', 'total');
            set(gcf1, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
            set(gca, 'FontSize', 16)
            set(gcf1, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            hold off;
        savefig(gcf1, figureName)

        bestThresholdTimeArr(q) = bestThresholdBreak;
        bestSwitchTimepointArr(q) = bestSwitchTime;

        figureName2 = sprintf('%sindividualMaxMaxB/concentration_x2x3ratio%s_gammaEqual.fig', folderName, num2str(ratio_2naught_3naught(q)));
        
        gcf2 = figure('Visible', 'off');
        hold on;
        plot(timeVector, concAVecBest, 'LineWidth', 3)
        plot(timeVector, concBVec, 'LineWidth', 3)
        set(gcf2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')')
        set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
        set(gca, 'FontSize', 16)
        hold off;
        savefig(gcf2, figureName2)
    end

    [failureTime, ind2] = max(bestThresholdTimeArr);
    switchpoint = bestSwitchTimepointArr(ind2);
end
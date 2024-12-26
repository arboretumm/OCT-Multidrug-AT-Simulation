function maxThreshBreak = multirun_uMaxA_uSingB_uneq_last(uMax, kR)
    
%population params that are staying the same

        %equation parameters
                %drug effect -- rn osimertinib and afatinib
                alphaA = 0.00155; %0.06; nM/hr
                alphaB = 0.00155; %0.06; nM/hr
            %population growth rate
                gammaN = 0.031; %0.031; %only EGFR+
                gammaA = 0.031;
                gammaB = 0.0155;
                gammaD = 0.0155; %0.011; %EGFR/T790M/C797S
                %no current difference between EGFR+ from L858R and ex19del
            %carrying capacity + population parameters
                kappa = 40000; %actual carrying capacity of cell culture // 40k 96MW
                thresholdRatio = kR;
                kappa_threshold = thresholdRatio.*kappa; %3500000; %goal capacity
            %drug limits
                uMaxA = uMax; %nM
                uMaxB = uMax; %nM
            %time parameters
                cellTime = 3600; %number of timesteps eqns are running for (in hours?)
            %initial conditions
                naivePopIC = 0.685.*kappa_threshold; 
                jointMutPopIC = 0.3.*kappa_threshold; 
                doubleMutPopIC = 0.015.*kappa_threshold;
            
    %ratios for looping through
        %ratio_2naught_3naught = [.05, .25, .5, .75, .95];
        ratio_2naught_3naught = [0.5];

        finalTimepoints = [];

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

            

            %feed new concentration vectors into the equations, then compare
            %previous simulation with new simulation & fail time
    
            initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];
        
            %solving/graphing system given parameters
            tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
            %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        
            Store.t=[];
            Store.x=[];
            t = tSPAN(1);
            x=initCond;
        
            if t(end) < tSPAN(end)
                options = odeset('Events', @(t,x) eventSubpopulationLimit(t, x, p, 5, (p.kappa_threshold+1)),'RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                [t, x] = ode45(@(t, x) controlAMax_controlBCalcThreshold(t, x, p), [t(end) tSPAN(end)], x(end,:), options);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
        
                if t(end) < tSPAN(end)
                    options3 = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4]);
                    [t, x] = ode45(@(t, x) controlAMax_controlBMax(t, x, p), [t(end) tSPAN(end)], x(end,:), options3);
                    Store.t=[Store.t ;t(2:end)];
                    Store.x=[Store.x; x(2:end,:)];
                end
            end

            figure(q)
            hold on
            plot(Store.t, Store.x, 'LineWidth', 3)
            linestyleorder(["-", '-', ':', '-', '-'])
            set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
            set(gca, 'FontSize', 16)
            hold off

            %find breach by pop 4 (first threshold break version)
            [rowsz, colsz] = size(Store.x);
            foundBreachC = false;
            foundBreachD = false;
            controlFailureTrack1 = 0;
            controlFailureTrack2 = 0;
            errorBreak = 1.001*(p.kappa_threshold);
            
            cf1Index = find((Store.x(:, 4) <= errorBreak), 1, 'last');
            cf2Index = find((Store.x(:, 5) <= errorBreak), 1, 'last');

            controlFailureTrack1 = Store.t(cf1Index);
            controlFailureTrack2 = Store.t(cf2Index);
    
            thresholdBreak = controlFailureTrack2;
            finalTimepoints(q) = thresholdBreak;
    end

    maxThreshBreak = max(finalTimepoints, [],'all');
end
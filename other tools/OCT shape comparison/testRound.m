function testRound(uMax, kR)

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
                naivePopIC = kappa_threshold - 10000;
            
    %ratios for looping through
        ratio_2naught_3naught = [.05, .25, .5, .75, .95];
        %ratio_2naught_3naught = [0.5]

        initCond = [naivePopIC 5000 5000 0 kappa_threshold];

        p = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, gammaN, gammaA, gammaB, gammaD, kappa, kappa_threshold, cellTime, naivePopIC, 5000, 5000, 0, uMaxA, uMaxB);
            % save(populationFileName, "popStruct");

        %solving/graphing system given parameters
            tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
            %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
        
            Store.t=[];
            Store.x=[];
            t = tSPAN(1);
            x=initCond;
        
            options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
            [t, x] = ode45(@(t, x) controlACalcThreshold_controlBMax(t, x, p), tSPAN, x(end,:), options);
            Store.t=[Store.t ;t(2:end)];
            Store.x=[Store.x; x(2:end,:)];

            figure(1)
            hold on
            plot(Store.t, Store.x, 'LineWidth', 3)
            linestyleorder(["-", '-', ':', '-', '-'])
            set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
            set(gca, 'FontSize', 16)
            hold off
end
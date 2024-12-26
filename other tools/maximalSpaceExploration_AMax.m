function maximalSpaceExploration_AMax()

%9/11/24 FOR LATER
%need to set it up so that it also goes through different x2:x3 ratios, try
%graphing the x2:x3 as separate graphs and as combined??

uMaxArr = [22.5, 45, 75, 150, 400; 15, 30, 50, 100, 400; 6, 12, 20, 40, 400];
[rowuMax, coluMax] = size(uMaxArr);
threshVec = [0.25, 0.5, 0.8];

%ratios for looping through
    ratio_2naught_3naught = [.05, .25, .5, .75, .95];
    %starting with just 0.5 though, change 3 below back to j for a loop

endTime = 1000 + 1;
timestep = 10;
timeVector = 1:timestep:endTime;

for q = 1:length(ratio_2naught_3naught)
    legendNames = [];
    failureArray = [];
    for k = 1:rowuMax
        for j = 1:coluMax

            kR = threshVec(k);
            uMax = uMaxArr(k, j);

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
                cellTime = 1000; %number of timesteps eqns are running for (in hours?)
            %initial conditions
                naivePopIC = 0.685.*kappa_threshold; 
                jointMutPopIC = 0.3.*kappa_threshold; 
                doubleMutPopIC = 0.015.*kappa_threshold;
        
        
        %will start a loop for ratio here i think
            
            %establish population
            initCondA = (1 - ratio_2naught_3naught(q)).*jointMutPopIC;
            initCondB = ratio_2naught_3naught(q).*jointMutPopIC;
        
            p = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, gammaN, gammaA, gammaB, gammaD, kappa, kappa_threshold, cellTime, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
            
            
            concAVec = uMaxA.*ones(1, length(timeVector));
            concBVec = zeros(1, length(timeVector));
        
            failureVec = zeros(1, length(timeVector));
        
        %     for i=1:length(timeVector)
            for i=1:length(timeVector)
                concBVecAlt = concBVec;
                concBVecAlt(i:end) = uMaxB;
        
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
                [t, x] = ode45(@(t, x) practicalityEquations(t, x, p, timeVector, concAVec, concBVecAlt), tSPAN, x(end,:), options);
                Store.t=[Store.t ;t(2:end)];
                Store.x=[Store.x; x(2:end,:)];
        
        %         figure(i)
        %         plot(Store.t, Store.x)
        
        
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
            failureArray = [failureArray; failureVec];
            stringName = sprintf('Thresh - %s, uMax - %s', num2str(threshVec(k)), num2str(uMaxArr(k, j)));
            legendNames = [legendNames, {stringName}];
        end
    end
    
%     figure(q)
%     hold on;
%     newcolorArr = ["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"];
%     plot(timeVector, failureArray(1:5, :), 'LineWidth', 2);
%     colororder(newcolorArr);
%     plot(timeVector, failureArray(6:10, :), 'LineWidth', 2, 'LineStyle','--');
%     colororder(newcolorArr);
%     plot(timeVector, failureArray(11:15, :), 'LineWidth', 2, 'LineStyle','-.');
%     colororder(newcolorArr);
%     hold off;
%     xlabel('time at switch')
%     ylabel('time at failure (first)')
%     title(sprintf("x2:x3 ratio - %s", num2str(ratio_2naught_3naught(q))))
%     legend(legendNames)
%     axis([0 1000 0 inf])
%     
%     set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
%     set(gca, 'FontSize', 16)

    figure(q)
    hold on;
    newcolorArr = ["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"];
    newLinePlan = ["-","--","-.",":", "-d"];
    plot(timeVector, failureArray(1:5, :), 'LineWidth', 2, 'Color', "#0072BD", 'MarkerIndices',(1:10:length(timeVector)));
    linestyleorder(newLinePlan);
    plot(timeVector, failureArray(6:10, :), 'LineWidth', 2, 'Color', "#EDB120", 'MarkerIndices', (1:10:length(timeVector)));
    linestyleorder(newLinePlan);
    plot(timeVector, failureArray(11:15, :), 'LineWidth', 2, 'Color', "#77AC30", 'MarkerIndices', (1:10:length(timeVector)));
    linestyleorder(newLinePlan);
    hold off;
    xlabel('time at switch')
    ylabel('time at failure (first)')
    title(sprintf("x2:x3 ratio - %s", num2str(ratio_2naught_3naught(q))))
    legend(legendNames)
    axis([0 1000 0 inf])
    
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
end
end
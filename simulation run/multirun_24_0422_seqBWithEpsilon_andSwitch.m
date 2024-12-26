function multirun_24_0422_seqBWithEpsilon_andSwitch(uMax, kR)
tic
%population set up
    %drug effect -- rn osimertinib and afatinib
        alphaA = 0.06; %0.06; nM/hr
        alphaB = 0.06; %0.06; nM/hr
    %population growth rate
        lambda_n = 0.031; %0.031; %only EGFR+
        baseJointLambda = 0.028; %0.022; %EGFR+/C797S
        lambda_d = 0.011; %0.011; %EGFR/T790M/C797S
        %no current difference between EGFR+ from L858R and ex19del
    %carrying capacity + population parameters
        kappa = 40000; %actual carrying capacity of cell culture // 40k 96MW
        thresholdRatio = kR;
        kappa_threshold = thresholdRatio.*kappa; %3500000; %goal capacity
    %drug limits
        uMaxA = uMax; %nM
        uMaxB = uMax; %nM
    %time parameters
        cellTime = 3000; %number of timesteps eqns are running for (in hours?)
    %initial conditions
        naivePopIC = 0.685.*kappa_threshold; 
        jointMutPopIC = 0.3.*kappa_threshold; 
        doubleMutPopIC = 0.015.*kappa_threshold;
    
%multirun parameters
    epsilon = 0.0075;
    ratio_2naught_3naught = [0.1, 0.2, 0.5, 1, 2, 5, 10];
    ratio_growthRates = [0.5, 0.75, 0.9, 1, 1.1111, 1.33333, 2];

    peakHeatmap_4 = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));
    peakHeatmap_5 = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));

folderName = sprintf("/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double boundary shape exploration/2024_0422_sequentialSwitch_Bfirst_uMax_thresh_tests_eps75/umax_%s_thresh_%s/", num2str(uMax), num2str(kR));

for i=1:length(ratio_2naught_3naught)
    for j=1:length(ratio_growthRates)
        %establish population
        popHold = jointMutPopIC./(1+(ratio_2naught_3naught(i)));
        initCondA = popHold;
        initCondB =(ratio_2naught_3naught(i)).*popHold;

        %establish growthRate
        gammaA = baseJointLambda;
        gammaB = ratio_growthRates(j).*baseJointLambda;

        %establish max - epsilon
        altMaxA = uMaxA - epsilon;
        altMaxB = uMaxB - epsilon;

        fileName = sprintf("individualGraphs/test_x2x3ratio_%s_gammaAgammaBratio_%s", num2str(ratio_2naught_3naught(i)), num2str(ratio_growthRates(j))); 

        popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, lambda_n, gammaA, gammaB, lambda_d, kappa, kappa_threshold, 3000, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
        [failureTime4, failureTime5] = sequentialDoubleBoundaryRunB_concDrugSwitch_withAOption(folderName, fileName, popStruct, altMaxA, altMaxB);
        idealTime = calculateIdealFromPopStruct(popStruct);
        holdTime4 = failureTime4/idealTime;
        holdTime5 = failureTime5/idealTime;
        %rows/y axis - growth rate, columns/x axis - population ratios
        peakHeatmap_4(j, i) = holdTime4;
        peakHeatmap_5(j, i) = holdTime5;
    end
end

    %graph heatmap
    figure('Visible', 'off');
    figureName = sprintf("Spread of Best (growth vs. pop ratio) for x(4)");
    imagesc(ratio_2naught_3naught, ratio_growthRates, peakHeatmap_4)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    %set(c3, 'ylim', [0 1]);
    caxis([0 1])
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    set(gcf, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    fileNamePeak2 = "heatmap_populationSpread_growth_x4";
    strFigSave5 = sprintf('%s%s.png', folderName, fileNamePeak2);
    strFigSave6 = sprintf('%s%s.fig', folderName, fileNamePeak2);
    saveas(gcf,strFigSave5); %saves population graph as fig
    saveas(gcf,strFigSave6); %saves population graph as png

    %graph heatmap
    figure('Visible', 'off');
    figureName = sprintf("Spread of Best (growth vs. pop ratio) for x(5)");
    imagesc(ratio_2naught_3naught, ratio_growthRates, peakHeatmap_5)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    %set(c3, 'ylim', [0 1]);
    caxis([0 1])
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    set(gcf, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    fileNamePeak2 = "heatmap_populationSpread_growth_x5";
    strFigSave5_2 = sprintf('%s%s.png', folderName, fileNamePeak2);
    strFigSave6_2 = sprintf('%s%s.fig', folderName, fileNamePeak2);
    saveas(gcf,strFigSave5_2); %saves population graph as fig
    saveas(gcf,strFigSave6_2); %saves population graph as png

toc
end
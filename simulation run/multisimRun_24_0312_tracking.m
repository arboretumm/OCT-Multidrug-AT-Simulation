function multisimRun_24_0312_tracking()
%FIGURE 1 FROM ELSA BOTH CONTROLS EMAIL
%goal: compare ratio of growth rates vs. ratio of concentrations @
%different values of x(2)naught/x(3)naught

%population params that are staying the same

%equation parameters
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
            kappa_threshold = 32000; %3500000; %goal capacity

    %drug limits
            uMaxA = 10000; %nM
            uMaxB = 10000; %nM

    %time parameters
        cellTime = 3000; %number of timesteps eqns are running for (in hours?)

    %initial conditions
        naivePopIC = 21500;
        jointMutPopIC = 10000;
        doubleMutPopIC = 500;


%ratios for looping through
ratio_2naught_3naught = [0.1, 0.2, 0.5, 1, 2, 5, 10];
ratio_concA_concB = [0.1, 0.2, 0.5, 1, 2, 5, 10]; %Cb/Ca
ratio_growthRates = [0.5, 0.75, 0.9, 1, 1.1111, 1.33333, 2];
peakHeatmap = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));

folderName = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double boundary shape exploration/2024_0312_populationRatioTests/";

%loop through simulations 
%reminder, for heatmaps - columns of array = x, rows of array = y
for i=1:length(ratio_2naught_3naught)
    heatmap = zeros(length(ratio_growthRates), length(ratio_concA_concB));
    fileNamePop = sprintf("test_x2x3ratio_%d_heatmap", ratio_2naught_3naught(i));
    for j=1:length(ratio_growthRates)
        for k=1:length(ratio_concA_concB)
            %establish population
            popHold = jointMutPopIC./(1+(ratio_2naught_3naught(i)));
            initCondA = popHold;
            initCondB =(ratio_2naught_3naught(i)).*popHold;

            %establish growthRate
            gammaA = baseJointLambda;
            gammaB = ratio_growthRates(j).*baseJointLambda;

            fileName = sprintf("test_x2x3ratio_%d_gammaAgammaBratio_%d_concAconcBratio_%d", ratio_2naught_3naught(i), ratio_growthRates(j), ratio_concA_concB(k)); 

            popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, lambda_n, gammaA, gammaB, lambda_d, kappa, kappa_threshold, 3000, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
            failureTime = jointDoubleBoundaryRun_constantProportionControls(folderName,fileName,popStruct, ratio_concA_concB(k));
            idealTime = calculateIdealFromPopStruct(popStruct);
            holdTime = failureTime/idealTime;
            heatmap(j,k) = holdTime;
            refTime = peakHeatmap(j, i);
            if refTime >= holdTime
                peakHeatmap(j, i) = refTime;
            else
                peakHeatmap(j, i) = holdTime;
            end
        end
    end
    %graph heatmap
    figureName = sprintf("Ratio of Simulated Failure and Ideal Failure - Ratio of x(2) and x(3) %d", ratio_2naught_3naught(i));
    imagesc(ratio_growthRates, ratio_concA_concB, heatmap)
    xlabel("\gamma_A : \gamma_B", 'FontSize', 18)
    ylabel("C_B/C_A", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    xtickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    ytickSpace = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
    xticks(xtickSpace);
    yticks(ytickSpace);
    xticklabels(ratio_growthRates);
    yticklabels(ratio_concA_concB);
    c3 = colorbar;
    c3.Label.String = "Ratio";
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)

    strFigSave1 = sprintf('%s%s.png', folderName, fileNamePop);
    strFigSave2 = sprintf('%s%s.fig', folderName, fileNamePop);
    saveas(gcf,strFigSave1); %saves population graph as fig
    saveas(gcf,strFigSave2); %saves population graph as png
end

%graph heatmap
    figureName = sprintf("Spread of Best");
    imagesc(ratio_growthRates, ratio_2naught_3naught, peakHeatmap)
    xlabel("\gamma_A : \gamma_B", 'FontSize', 18)
    ylabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    c3.Label.String = "Ratio";
    xtickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    ytickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    xticks(xtickSpace);
    yticks(ytickSpace);
    xticklabels(ratio_growthRates);
    yticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    fileNamePeak = "heatmap_populationSpread";
    strFigSave3 = sprintf('%s%s.png', folderName, fileNamePeak);
    strFigSave4 = sprintf('%s%s.fig', folderName, fileNamePeak);
    saveas(gcf,strFigSave3); %saves population graph as fig
    saveas(gcf,strFigSave4); %saves population graph as png

end
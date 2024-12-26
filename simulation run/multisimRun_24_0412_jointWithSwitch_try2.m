function multisimRun_24_0412_jointWithSwitch_try2(uMax, kR)
tic
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


%ratios for looping through
ratio_2naught_3naught = [0.1, 0.2, 0.5, 1, 2, 5, 10];
ratio_concA_concB = [0.1, 0.2, 0.5, 1, 2, 5, 10]; %Cb/Ca
ratio_growthRates = [0.5, 0.75, 0.9, 1, 1.1111, 1.33333, 2];
peakHeatmap_4 = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));
peakHeatmap_5 = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));
peakHeatmap2_4 = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));
peakHeatmap2_5 = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));

folderName = sprintf("/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double boundary shape exploration/2024_0412_jointWithSwitch2_uMax_thresh_tests/umax_%s_thresh_%s/", num2str(uMax), num2str(kR));

%loop through simulations 
%reminder, for heatmaps - columns of array = x, rows of array = y
for i=1:length(ratio_2naught_3naught)
    heatmap4 = zeros(length(ratio_growthRates), length(ratio_concA_concB));
    heatmap5 = zeros(length(ratio_growthRates), length(ratio_concA_concB));
    fileNamePop = sprintf("test_x2x3ratio_%s_heatmap", num2str(ratio_2naught_3naught(i)));
    for j=1:length(ratio_growthRates)
        for k=1:length(ratio_concA_concB)
            %establish population
            popHold = jointMutPopIC./(1+(ratio_2naught_3naught(i)));
            initCondA = popHold;
            initCondB =(ratio_2naught_3naught(i)).*popHold;

            %establish growthRate
            gammaA = baseJointLambda;
            gammaB = ratio_growthRates(j).*baseJointLambda;

            fileName = sprintf("individualGraphs/test_x2x3ratio_%s_gammaAgammaBratio_%s_concAconcBratio_%s", num2str(ratio_2naught_3naught(i)), num2str(ratio_growthRates(j)), num2str(ratio_concA_concB(k))); 

            popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, lambda_n, gammaA, gammaB, lambda_d, kappa, kappa_threshold, 3000, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
            [failureTime4, failureTime5] = jointDoubleBoundaryRun_constantProportionControlswithSwitchtry2(folderName,fileName,popStruct, ratio_concA_concB(k));
            idealTime = calculateIdealFromPopStruct(popStruct);
            holdTime4 = failureTime4/idealTime;
            holdTime5 = failureTime5/idealTime;
            heatmap4(j,k) = holdTime4;
            heatmap5(j,k) = holdTime5;
            refTime4 = peakHeatmap_4(j, i);
            refTime5 = peakHeatmap_5(j, i);
            if refTime4 >= holdTime4
                peakHeatmap_4(j, i) = refTime4;
            else
                peakHeatmap_4(j, i) = holdTime4;
            end
            if refTime5 >= holdTime5
                peakHeatmap_5(j, i) = refTime5;
            else
                peakHeatmap_5(j, i) = holdTime5;
            end
            refTime2_4 = peakHeatmap2_4(k, i);
            refTime2_5 = peakHeatmap2_5(k, i);
            if refTime2_4 >= holdTime4
                peakHeatmap2_4(k, i) = refTime2_4;
            else
                peakHeatmap2_4(k, i) = holdTime4;
            end
            if refTime2_5 >= holdTime5
                peakHeatmap2_5(k, i) = refTime2_5;
            else
                peakHeatmap2_5(k, i) = holdTime5;
            end
        end
    end
    %graph heatmap
    figureName = sprintf("Ratio of Simulated Failure and Ideal Failure for x(4) - Ratio of x(2) and x(3) %s", num2str(ratio_2naught_3naught(i)));
    imagesc(ratio_concA_concB, ratio_growthRates, heatmap4)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("C_B/C_A", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_concA_concB);
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)

    strFigSave1 = sprintf('%s%s_x4.png', folderName, fileNamePop);
    strFigSave2 = sprintf('%s%s_x4.fig', folderName, fileNamePop);
    saveas(gcf,strFigSave1); %saves population graph as fig
    saveas(gcf,strFigSave2); %saves population graph as png

    %graph heatmap
    figureName = sprintf("Ratio of Simulated Failure and Ideal Failure for x(5) - Ratio of x(2) and x(3) %s", num2str(ratio_2naught_3naught(i)));
    imagesc(ratio_concA_concB, ratio_growthRates, heatmap5)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("C_B/C_A", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_concA_concB);
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)

    strFigSave1 = sprintf('%s%s_x5.png', folderName, fileNamePop);
    strFigSave2 = sprintf('%s%s_x5.fig', folderName, fileNamePop);
    saveas(gcf,strFigSave1); %saves population graph as fig
    saveas(gcf,strFigSave2); %saves population graph as png

end

%graph heatmap
    figureName = sprintf("Spread of Best (concentration vs. pop ratio) for x(4)");
    imagesc(ratio_2naught_3naught, ratio_concA_concB, peakHeatmap2_4)
    ylabel("conc A : conc B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_concA_concB);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    fileNamePeak = "heatmap_populationSpread_concentration_x4";
    strFigSave3 = sprintf('%s%s.png', folderName, fileNamePeak);
    strFigSave4 = sprintf('%s%s.fig', folderName, fileNamePeak);
    saveas(gcf,strFigSave3); %saves population graph as fig
    saveas(gcf,strFigSave4); %saves population graph as png

    %graph heatmap
    figureName = sprintf("Spread of Best (concentration vs. pop ratio) for x(5)");
    imagesc(ratio_2naught_3naught, ratio_concA_concB, peakHeatmap2_5)
    ylabel("conc A : conc B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_concA_concB);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    fileNamePeak = "heatmap_populationSpread_concentration_x5";
    strFigSave3_2 = sprintf('%s%s.png', folderName, fileNamePeak);
    strFigSave4_2 = sprintf('%s%s.fig', folderName, fileNamePeak);
    saveas(gcf,strFigSave3_2); %saves population graph as fig
    saveas(gcf,strFigSave4_2); %saves population graph as png


  %graph heatmap
    figureName = sprintf("Spread of Best (growth vs. pop ratio) for x(4)");
    imagesc(ratio_2naught_3naught, ratio_growthRates, peakHeatmap_4)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    fileNamePeak2 = "heatmap_populationSpread_growth_x4";
    strFigSave5 = sprintf('%s%s.png', folderName, fileNamePeak2);
    strFigSave6 = sprintf('%s%s.fig', folderName, fileNamePeak2);
    saveas(gcf,strFigSave5); %saves population graph as fig
    saveas(gcf,strFigSave6); %saves population graph as png

     %graph heatmap
    figureName = sprintf("Spread of Best (growth vs. pop ratio) for x(5)");
    imagesc(ratio_2naught_3naught, ratio_growthRates, peakHeatmap_5)
    ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
    xlabel("x2:x3", 'FontSize', 18)
    axis xy;
    title(figureName, 'FontSize', 24)
    c3 = colorbar;
    caxis([0 1]);
    c3.Label.String = "Ratio";
    ytickSpace = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
    xtickSpace = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
    yticks(ytickSpace);
    xticks(xtickSpace);
    yticklabels(ratio_growthRates);
    xticklabels(ratio_2naught_3naught);
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    fileNamePeak2 = "heatmap_populationSpread_growth_x5";
    strFigSave5_2 = sprintf('%s%s.png', folderName, fileNamePeak2);
    strFigSave6_2 = sprintf('%s%s.fig', folderName, fileNamePeak2);
    saveas(gcf,strFigSave5_2); %saves population graph as fig
    saveas(gcf,strFigSave6_2); %saves population graph as png
toc
end

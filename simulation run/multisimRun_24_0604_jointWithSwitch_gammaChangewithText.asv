function multisimRun_24_0604_jointWithSwitch_gammaChangewithText(uMax, kR)
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

    %establish heatmaps -- all x(5) this time around 
        heatmap_gammaPop = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));
        heatmap_gammaPop_worstVal = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));
        heatmap_gammaPop_concVal = zeros(length(ratio_growthRates), length(ratio_2naught_3naught));

        heatmap_concPop = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));
        heatmap_concPop_worstVal = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));
        heatmap_concPop_gammaVal = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));
    
    folderName = sprintf("/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double boundary shape exploration/2024_0605_jointSwitch_withText/umax_%s_thresh_%s/", num2str(uMax), num2str(kR)); 

    %run heatmap loops
    for i=1:length(ratio_2naught_3naught)  
        heatmapVals = zeros(length(ratio_growthRates), length(ratio_concA_concB));
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
                normedTime = failureTime5./idealTime;

                %heatmaps
                heatmapVals(j,k) = normedTime;
                if normedTime > heatmap_gammaPop(j, i) %best time for gamma stuff
                    heatmap_gammaPop(j, i) = normedTime;
                    heatmap_gammaPop_concVal(j, i) = ratio_concA_concB(k);
                end
                if normedTime < heatmap_gammaPop_worstVal(j, i) || heatmap_gammaPop_worstVal(j, i) == 0
                    heatmap_gammaPop_worstVal(j, i) = normedTime;
                end

                if normedTime > heatmap_concPop(k, i)
                    heatmap_concPop(k, i) = normedTime;
                    heatmap_concPop_gammaVal(k, i) = ratio_growthRates(j);
                end
                if normedTime < heatmap_concPop_worstVal(k, i) || heatmap_concPop_worstVal(k, i) == 0
                    heatmap_concPop_worstVal(k, i) = normedTime;
                end
            end
        end

    %graph heatmap
    figureName = sprintf("Ratio of Simulated Failure and Ideal Failure for x(5) - Ratio of x(2) and x(3) %s", num2str(ratio_2naught_3naught(i)));
    imagesc(ratio_concA_concB, ratio_growthRates, heatmapVals)
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

    heatmap_gammaPop_concVal
    heatmap_gammaPop_worstVal
    heatmap_concPop_gammaVal
    heatmap_concPop_worstVal

    %out here, need to code the heatmaps blank & heatmaps with labels 

    %GAMMA V POP SECTION
        %blank heatmap gamma v pop
        figureName = sprintf("Spread of Best (growth vs. pop ratio)");
        imagesc(ratio_2naught_3naught, ratio_growthRates, heatmap_gammaPop)
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
        fileNamePeakgamma = "heatmap_populationSpread_growth_blank";
        strFigSavegamma_1a = sprintf('%s%s.png', folderName, fileNamePeakgamma);
        strFigSavegamma_1b = sprintf('%s%s.fig', folderName, fileNamePeakgamma);
        saveas(gcf,strFigSavegamma_1a); %saves population graph as fig
        saveas(gcf,strFigSavegamma_1b); %saves population graph as png

        
        rowArr = linspace(ratio_growthRates(1), ratio_growthRates(end), length(ratio_growthRates));
        colArr = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
        %heatmap gamma v pop with worst time labeled ALT OPT
        figureName = sprintf("Spread of Best (growth vs. pop ratio), worst time");
        hold on
        imagesc(ratio_2naught_3naught, ratio_growthRates, heatmap_gammaPop)
        for i=1:length(ratio_2naught_3naught)
            for j=1:length(ratio_growthRates)
                num2str(heatmap_gammaPop_worstVal(j, i))
                text(colArr(i), rowArr(j), num2str(heatmap_gammaPop_worstVal(j, i)));
            end
        end
        hold off
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
        fileNamePeakgamma_worst = "heatmap_populationSpread_growth_worstLabel";
        strFigSavegamma_3a = sprintf('%s%s.png', folderName, fileNamePeakgamma_worst);
        strFigSavegamma_3b = sprintf('%s%s.fig', folderName, fileNamePeakgamma_worst);
        saveas(gcf,strFigSavegamma_3a); %saves population graph as fig
        saveas(gcf,strFigSavegamma_3b); %saves population graph as png

        %heatmap gamma v pop with chosen conc labeled 
        figureName = sprintf("Spread of Best (growth vs. pop ratio), conc label");
        hold on
        imagesc(ratio_2naught_3naught, ratio_growthRates, heatmap_gammaPop)
        for i=1:length(ratio_2naught_3naught)
            for j=1:length(ratio_growthRates)
                text(colArr(i), rowArr(j), num2str(heatmap_gammaPop_concVal(j, i)));
            end
        end
        hold off
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
        fileNamePeakgamma_conc = "heatmap_populationSpread_growth_concLabel";
        strFigSavegamma_4a = sprintf('%s%s.png', folderName, fileNamePeakgamma_conc);
        strFigSavegamma_4b = sprintf('%s%s.fig', folderName, fileNamePeakgamma_conc);
        saveas(gcf,strFigSavegamma_4a); %saves population graph as fig
        saveas(gcf,strFigSavegamma_4b); %saves population graph as png

    %CONC V POP SECTION

        %blank heatmap conc v pop
        figureName = sprintf("Spread of Best (conc vs. pop ratio)");
        imagesc(ratio_2naught_3naught, ratio_concA_concB, heatmap_concPop)
        ylabel("C_B:C_A", 'FontSize', 18)
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
        fileNamePeakConc = "heatmap_populationSpread_conc_blank";
        strFigSaveConc_1a = sprintf('%s%s.png', folderName, fileNamePeakConc);
        strFigSaveConc_1b = sprintf('%s%s.fig', folderName, fileNamePeakConc);
        saveas(gcf,strFigSaveConc_1a); %saves population graph as fig
        saveas(gcf,strFigSaveConc_1b); %saves population graph as png

        
        rowArr2 = linspace(ratio_concA_concB(1), ratio_concA_concB(end), length(ratio_concA_concB));
        colArr2 = linspace(ratio_2naught_3naught(1), ratio_2naught_3naught(end), length(ratio_2naught_3naught));
        %heatmap conc v pop with worst time labeled ALT OPT
        figureName = sprintf("Spread of Best (convc vs. pop ratio), worst time");
        hold on
        imagesc(ratio_2naught_3naught, ratio_concA_concB, heatmap_concPop)
        for i=1:length(ratio_2naught_3naught)
            for j=1:length(ratio_concA_concB)
                num2str(heatmap_concPop_worstVal(j, i))
                text(colArr2(i), rowArr2(j), num2str(heatmap_concPop_worstVal(j, i)));
            end
        end
        hold off
        ylabel("C_B:C_A", 'FontSize', 18)
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
        fileNamePeakConc_worst = "heatmap_populationSpread_conc_worstLabel";
        strFigSaveConc_2a = sprintf('%s%s.png', folderName, fileNamePeakConc_worst);
        strFigSaveConc_2b = sprintf('%s%s.fig', folderName, fileNamePeakConc_worst);
        saveas(gcf,strFigSaveConc_2a); %saves population graph as fig
        saveas(gcf,strFigSaveConc_2b); %saves population graph as png

        %heatmap gamma v pop with chosen conc labeled 
        figureName = sprintf("Spread of Best (conc vs. pop ratio), gamma");
        hold on
        imagesc(ratio_2naught_3naught, ratio_concA_concB, heatmap_concPop)
        for i=1:length(ratio_2naught_3naught)
            for j=1:length(ratio_concA_concB)
                num2str(heatmap_concPop_gammaVal(j, i))
                text(colArr2(i), rowArr2(j), num2str(heatmap_concPop_gammaVal(j, i)));
            end
        end
        hold off
        ylabel("\gamma_A : \gamma_B", 'FontSize', 18)
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
        fileNamePeakConc_gamma = "heatmap_populationSpread_conc_gammaLabel";
        strFigSaveConc_3a = sprintf('%s%s.png', folderName, fileNamePeakConc_gamma);
        strFigSaveConc_3b = sprintf('%s%s.fig', folderName, fileNamePeakConc_gamma);
        saveas(gcf,strFigSaveConc_3a); %saves population graph as fig
        saveas(gcf,strFigSaveConc_3b); %saves population graph as png

    toc
end

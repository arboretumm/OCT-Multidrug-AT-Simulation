function multisimRun_24_0828_jointWithSwitch_heatmapSave_GammaEqual_last(uMax, kR, folderBase)
tic
%FIGURE 1 FROM ELSA BOTH CONTROLS EMAIL
%goal: compare ratio of growth rates vs. ratio of concentrations @
%different values of x(2)naught/x(3)naught

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
        cellTime = 3000; %number of timesteps eqns are running for (in hours?)
    %initial conditions
        naivePopIC = 0.685.*kappa_threshold; 
        jointMutPopIC = 0.3.*kappa_threshold; 
        doubleMutPopIC = 0.015.*kappa_threshold;

    %ratios for looping through
        ratio_2naught_3naught = [.05, .25, .5, .75, .95];
        ratio_concA_concB = [0.1, 0.2, 1, 5, 10];
    
    folderName = sprintf("%sjoint equal/umax_%s_thresh_%s/", folderBase, num2str(uMax), num2str(kR));
    
    vectorName = sprintf("%sVectorRange.mat", folderName);
    save(vectorName, "ratio_2naught_3naught", "ratio_concA_concB");

    heatmapVals = zeros(length(ratio_concA_concB), length(ratio_2naught_3naught));

    for i=1:length(ratio_2naught_3naught)
        for j=1:length(ratio_concA_concB)
            %establish population
            initCondA = (1 - ratio_2naught_3naught(i)).*jointMutPopIC;
            initCondB = ratio_2naught_3naught(i).*jointMutPopIC;
            
            %growth rate already established

            %pop struct & file names
            fileName = sprintf("individualGraphs/test_x2x3ratio%s_concAconcBratio%s_gammaEqual", num2str(ratio_2naught_3naught(i)), num2str(ratio_concA_concB(j)));
            populationFileName = sprintf("%sindividualGraphs/population_x2x3ratio%s_concAconcBratio%s_gammaEqual.mat", folderName, num2str(ratio_2naught_3naught(i)), num2str(ratio_concA_concB(j)));
            popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, gammaN, gammaA, gammaB, gammaD, kappa, kappa_threshold, 1000, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
            save(populationFileName, "popStruct");

            %results
            [failureTime4, failureTime5] = jointDoubleBoundaryRun_constantProportionControls_last(folderName,fileName,popStruct, ratio_concA_concB(j));
            idealTime = calculateIdealFromPopStruct(popStruct);
            %normedTime = failureTime5./idealTime;

            %heatmaps
            %heatmapVals(j,i) = normedTime; %rows conc and columns pop?
            heatmapVals(j, i) = failureTime5;

        end
    end
    heatmapFile = sprintf("%ssimHeatmap.mat", folderName);
    save(heatmapFile, "heatmapVals");

    toc
end

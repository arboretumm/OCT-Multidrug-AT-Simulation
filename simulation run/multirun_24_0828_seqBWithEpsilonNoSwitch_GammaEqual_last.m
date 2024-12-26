function multirun_24_0828_seqBWithEpsilonNoSwitch_GammaEqual_last(uMax, kR, folderBase)
tic
%population set up
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
    
%multirun parameters
    epsilon = 0.000;
    ratio_2naught_3naught = [.05, .25, .5, .75, .95];


folderName = sprintf("%ssequential equal/umax_%s_thresh_%s/", folderBase, num2str(uMax), num2str(kR));

heatvectorVals = zeros(length(ratio_2naught_3naught));

vectorName = sprintf("%sVectorRange.mat", folderName);
save(vectorName, "ratio_2naught_3naught", "epsilon");

    for i=1:length(ratio_2naught_3naught)
        %establish population
        initCondA = (1 - ratio_2naught_3naught(i)).*jointMutPopIC;
        initCondB = ratio_2naught_3naught(i).*jointMutPopIC;
        
        %growth rate already established
    
        %uMax alteration
        altMaxA = uMaxA - epsilon;
        altMaxB = uMaxB - epsilon;
    
        %pop struct & file names
        fileName = sprintf("individualGraphs/test_x2x3ratio%s_gammaEqual", num2str(ratio_2naught_3naught(i)));
        populationFileName = sprintf("%sindividualGraphs/population_x2x3ratio%s_gammaEqual.mat", folderName, num2str(ratio_2naught_3naught(i)));
        popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, gammaN, gammaA, gammaB, gammaD, kappa, kappa_threshold, 3000, naivePopIC, initCondA, initCondB, doubleMutPopIC, altMaxA, altMaxB);
        save(populationFileName, "popStruct");
    

        [failureTime4, failureTime5] = sequentialDoubleBoundaryRunB_concDrugSwitch_Aopt_last(folderName, fileName, popStruct, altMaxA, altMaxB);
        idealTime = calculateIdealFromPopStruct(popStruct);
        %normedTime = failureTime5./idealTime;
    
        %heatvectorVals(i) = normedTime;
        heatvectorVals(i) = failureTime5;
    end

    heatmapFile = sprintf("%ssimHeatvector.mat", folderName);
    save(heatmapFile, "heatvectorVals");

toc
end
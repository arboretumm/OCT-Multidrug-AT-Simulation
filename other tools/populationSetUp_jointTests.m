function popStruct = populationSetUp_jointTests(popRatio, gammaRatio, kR, uMax)
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

    %popRatio
        popHold = jointMutPopIC./(1+popRatio);
        initCondA = popHold;
        initCondB = popRatio.*popHold;

    %gammaRatio
        gammaA = baseJointLambda;
        gammaB = gammaRatio.*baseJointLambda;

    %population
    popStruct = fullEditPopulationParametersWithDrugLimits(alphaA, alphaB, lambda_n, gammaA, gammaB, lambda_d, kappa, kappa_threshold, 3000, naivePopIC, initCondA, initCondB, doubleMutPopIC, uMaxA, uMaxB);
end

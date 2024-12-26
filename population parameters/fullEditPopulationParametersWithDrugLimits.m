function s = fullEditPopulationParametersWithDrugLimits(aA, aB, lN, lM, lS, lD, k, kT, time, nPop, aPop, bPop, dPop, cMaxA, cMaxB)
    %DESCRIPTION OF PARAMETER BUILD
        %no drug-based suppression/competition
        %bondarenko + chmielecki derived terms in commented section
        %OTHER?

%parameter definition
    %equation parameters
        %drug effect -- rn osimertinib and afatinib
            s.alphaA = aA; %0.06; nM/hr
            s.alphaB = aB; %0.06; nM/hr
        %population growth rate
            s.lambda_n = lN; %0.031; %only EGFR+
            s.lambda_m = lM; %0.022; %EGFR+/T790M
            s.lambda_s = lS;%0.022; %EGFR+/C797S
            s.lambda_d = lD;%0.011; %EGFR/T790M/C797S
            %no current difference between EGFR+ from L858R and ex19del
        %carrying capacity + population parameters
            s.kappa = k; %4800000; %actual carrying capacity of cell culture // 40k 96MW
            s.kappa_threshold = kT; %3500000; %goal capacity

    %drug limits
            s.uMaxA = cMaxA; %nM ^(-1)
            s.uMaxB = cMaxB; %nM ^(-1)

    %time parameters
        s.cellTime = time; %number of timesteps eqns are running for (in hours?)
        s.timeRatio = 1000; %how many divisions each timestep will be split into
        s.timeChunk = s.cellTime.*s.timeRatio; %individual chunks of time used later

    %initial conditions
        s.naivePopIC = nPop;
        s.mutAPopIC = aPop;
        s.mutBPopIC = bPop;
        s.doubleMutPopIC = dPop;
        s.totalPopIC = s.mutAPopIC + s.mutBPopIC + s.doubleMutPopIC + s.naivePopIC;

end
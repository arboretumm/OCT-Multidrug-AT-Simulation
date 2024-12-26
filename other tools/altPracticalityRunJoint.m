function altPracticalityRunJoint(popStructFile, practicalSave, resultsFile, constProp, timestep, dosingInc)
    load(popStructFile, "popStruct"); 
    p = popStruct;

    initCond = [p.naivePopIC, p.mutAPopIC, p.mutBPopIC, p.doubleMutPopIC, p.totalPopIC];
    
    timeChunks = floor(p.cellTime./timestep);
    endTime = timeChunks*timestep;

    timeVals = 0:timestep:endTime;
    concentrationARecord = zeros(1, length(timeVals));
    concentrationBRecord = zeros(1, length(timeVals));

    %solving/graphing system given parameters
    tSPAN = [1 p.cellTime]; %timespan the ODE is being solved over
    %options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);

    Store.t=[];
    Store.x=[];
    x=initCond;
    
    concentrationAStart = ((1 - p.kappa_threshold/p.kappa)*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC))./(p.alphaA.*(p.naivePopIC+p.mutAPopIC)+p.alphaB.*constProp.*(p.naivePopIC+p.mutBPopIC))
    concentrationBStart = constProp*concentrationAStart

    if concentrationAStart < 0
        concentrationAStart = 0;
        recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC) - (p.alphaA.*concentrationAStart.*(p.naivePopIC+p.mutAPopIC)))./(p.alphaB.*(p.naivePopIC+p.mutBPopIC)); %insertequationusingconcAinstead
        if recalcB < 0
            concentrationBStart = 0;
        elseif recalcB > p.uMaxB
            concentrationBStart = p.uMaxB;
        else
            concentrationBStart = recalcB;
        end
    elseif concentrationAStart > p.uMaxA
        concentrationAStart = p.uMaxA;
        concA = p.uMaxA;
        recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC) - (p.alphaA.*concA.*(p.naivePopIC+p.mutAPopIC)))./(p.alphaB.*(p.naivePopIC+p.mutBPopIC)); %insertequationusingconcAinstead
        if recalcB < 0
            concentrationBStart = 0;
        elseif recalcB > p.uMaxB
            concentrationBStart = p.uMaxB;
        else
            concentrationBStart = recalcB;
        end
    else
        if concentrationBStart < 0
            concentrationBStart = 0;
            concB = 0;
            recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC) - (p.alphaB.*concB.*(p.naivePopIC+p.mutBPopIC)))./(p.alphaA.*(p.naivePopIC+p.mutAPopIC));
            if recalcA < 0
                concentrationAStart = 0;
            elseif recalcA > p.uMaxA
                concentrationAStart = p.uMaxA;
            else
                concentrationAStart = recalcA;
            end
        elseif concentrationBStart > p.uMaxB
            concentrationBStart = p.uMaxB;
            concB = p.uMaxB;
            recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* p.naivePopIC + p.lambda_m .* p.mutAPopIC + p.lambda_s .* p.mutBPopIC + p.lambda_d .* p.doubleMutPopIC) - (p.alphaB.*concB.*(p.naivePopIC+p.mutBPopIC)))./(p.alphaA.*(p.naivePopIC+p.mutAPopIC));
            if recalcA < 0
                concentrationAStart = 0;
            elseif recalcA > p.uMaxA
                concentrationAStart = p.uMaxA;
            else
                concentrationAStart = recalcA;
            end
        else
            concentrationAStart = concentrationAStart;
            concentrationBStart = concentrationBStart;
        end
    end

    roundedConcAStart = ceil(concentrationAStart/dosingInc)*dosingInc;
    roundedConcBStart = ceil(concentrationBStart/dosingInc)*dosingInc;
    concentrationARecord(1) = roundedConcAStart;
    concentrationBRecord(1) = roundedConcBStart;

    options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3 4 5]);
    [t, x] = ode45(@(t, x) practicalityEquationsTimestep(t, x, p, roundedConcAStart, roundedConcBStart), [0 timestep], x(end,:), options);
    Store.t=[Store.t ;t(2:end)];
    Store.x=[Store.x; x(2:end,:)];
    
    n = 2;

    while Store.t < tSPAN(end)
        concA = ((1 - p.kappa_threshold/p.kappa)*(p.lambda_n .* x(end,1) + p.lambda_m .* x(end,2) + p.lambda_s .* x(end,3) + p.lambda_d .* x(end,4)))./(p.alphaA.*(x(end,1)+x(end,2))+p.alphaB.*constProp.*(x(end,1)+x(end,3)));
        concB = constProp*concA;
        if concA < 0
            concA = 0;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(end,1) + p.lambda_m .* x(end,2) + p.lambda_s .* x(end,3) + p.lambda_d .* x(end,4)) - (p.alphaA.*concA.*(x(end,1)+x(end,2))))./(p.alphaB.*(x(end, 1)+x(end, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concB = 0;
            elseif recalcB > p.uMaxB
                concB = p.uMaxB;
            else
                concB = recalcB;
            end
        elseif concA > p.uMaxA
            concA = p.uMaxA;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(end,1) + p.lambda_m .* x(end,2) + p.lambda_s .* x(end,3) + p.lambda_d .* x(end,4)) - (p.alphaA.*concA.*(x(end,1)+x(end,2))))./(p.alphaB.*(x(end, 1)+x(end, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concB = 0;
            elseif recalcB > p.uMaxB
                concB = p.uMaxB;
            else
                concB = recalcB;
            end
        else
            if concB < 0
                concB = 0;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(end,1) + p.lambda_m .* x(end,2) + p.lambda_s .* x(end,3) + p.lambda_d .* x(end,4)) - (p.alphaB.*concB.*(x(end,1)+x(end,3))))./(p.alphaA.*(x(end,1)+x(end,2)));
                if recalcA < 0
                    concA = 0;
                elseif recalcA > p.uMaxA
                    concA = p.uMaxA;
                else
                    concA = recalcA;
                end
            elseif concB > p.uMaxB
                concB = p.uMaxB;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* x(end,1) + p.lambda_m .* x(end,2) + p.lambda_s .* x(end,3) + p.lambda_d .* x(end,4)) - (p.alphaB.*concB.*(x(end,1)+x(end,3))))./(p.alphaA.*(x(end,1)+x(end,2)));
                if recalcA < 0
                    concA = 0;
                elseif recalcA > p.uMaxA
                    concA = p.uMaxA;
                else
                    concA = recalcA;
                end
            else
                concA = concA;
                concB = concB;
            end
        end
        roundConcentrationA = ceil(concA/dosingInc)*dosingInc;
        roundConcentrationB = ceil(concB/dosingInc)*dosingInc;
        concentrationARecord(n) = roundConcentrationA;
        concentrationBRecord(n) = roundConcentrationB;

        [t, x] = ode45(@(t,x) practicalityEquationsTimestep(t, x, p, roundConcentrationA, roundConcentrationB), [(n-1)*timestep, (n)*timestep], x(end, :), options);
        Store.t=[Store.t ;t(2:end)];
        Store.x=[Store.x; x(2:end,:)];

        n = n+1;
    end

    load(resultsFile, "tResult");
    load(resultsFile, "xResult");

    %back calculate the concentration of the original run -- 
        concentrationAOld = ((1 - p.kappa_threshold/p.kappa)*(p.lambda_n .* xResult(:,1) + p.lambda_m .* xResult(:,2) + p.lambda_s .* xResult(:,3) + p.lambda_d .* xResult(:,4)))./(p.alphaA.*(xResult(:,1)+xResult(:,2))+p.alphaB.*constProp.*(xResult(:,1)+xResult(:,3)));
        concentrationBOld = constProp*concentrationAOld;

        concentrationA = zeros(1, length(concentrationAOld));
        concentrationB = zeros(1, length(concentrationBOld));

    for i=1:length(concentrationAOld)
        concA = concentrationAOld(i);
        concB = concentrationBOld(i);
        if concA < 0
            concentrationA(1, i) = 0;
            concA = 0;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaA.*concA.*(xResult(i,1)+xResult(i,2))))./(p.alphaB.*(xResult(i, 1)+xResult(i, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concentrationB(1, i) = 0;
            elseif recalcB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
            else
                concentrationB(1, i) = recalcB;
            end
        elseif concA > p.uMaxA
            concentrationA(1, i) = p.uMaxA;
            concA = p.uMaxA;
            recalcB = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaA.*concA.*(xResult(i,1)+xResult(i,2))))./(p.alphaB.*(xResult(i, 1)+xResult(i, 3))); %insertequationusingconcAinstead
            if recalcB < 0
                concentrationB(1, i) = 0;
            elseif recalcB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
            else
                concentrationB(1, i) = recalcB;
            end
        else
            if concB < 0
                concentrationB(1, i) = 0;
                concB = 0;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaB.*concB.*(xResult(i,1)+xResult(i,3))))./(p.alphaA.*(xResult(i,1)+xResult(i,2)));
                if recalcA < 0
                    concentrationA(1, i) = 0;
                elseif recalcA > p.uMaxA
                    concentrationA(1, i) = p.uMaxA;
                else
                    concentrationA(1, i) = recalcA;
                end
            elseif concB > p.uMaxB
                concentrationB(1, i) = p.uMaxB;
                concB = p.uMaxB;
                recalcA = ((1 - p.kappa_threshold./p.kappa).*(p.lambda_n .* xResult(i,1) + p.lambda_m .* xResult(i,2) + p.lambda_s .* xResult(i,3) + p.lambda_d .* xResult(i,4)) - (p.alphaB.*concB.*(xResult(i,1)+xResult(i,3))))./(p.alphaA.*(xResult(i,1)+xResult(i,2)));
                if recalcA < 0
                    concentrationA(1, i) = 0;
                elseif recalcA > p.uMaxA
                    concentrationA(1, i) = p.uMaxA;
                else
                    concentrationA(1, i) = recalcA;
                end
            else
                concentrationA(1, i) = concA;
                concentrationB(1, i) = concB;
            end
        end
    end

    gcf1 = figure(2);
    hold on
    plot(timeVals, concentrationARecord, 'LineWidth', 3, 'LineStyle', '-', 'Color', '#0072BD')
    plot(tResult, concentrationA, 'LineWidth', 3, 'LineStyle', '--', 'Color', '#9fc8e2')


    plot(timeVals, concentrationBRecord, 'LineWidth', 3, 'LineStyle', '-', 'Color', '#c4a55c')
    plot(tResult, concentrationB, 'LineWidth', 3, 'LineStyle', '--', 'Color', '#EDB120')

    title('Plot of drug concentrations', 'FontSize', 24);
    xlabel('time (hours)', 'FontSize', 18);
    ylabel('concentration (nM)', 'FontSize', 18);
    legend('[drug A]', 'practical [drug A]', '[drug B]', 'practical [drug B]');
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    hold off
    %save figure
    strFileName1 = sprintf('%s%s.fig', practicalSave, "_drugConcentrationApprox");
    strFileName2 = sprintf('%s%s.png', practicalSave, "_drugConcentrationApprox");
    saveas(gcf1,strFileName1); %saves population graph as fig
    saveas(gcf1,strFileName2); %saves population graph as png
     

%     gcf2 = figure('Visible', 'off');
    gcf2 = figure(1);
    hold on;
    line_color = ['#0072BD', 	"#D95319", 	"#EDB120", 	"#7E2F8E", 	"#77AC30", "#9fc8e2", 	"#d8b7a3", 	"#c4a55c", 	"#c9b4cd", 	"#bfcdac"];
    plot(Store.t, Store.x, 'LineWidth', 2.5);
    plot(tResult, xResult, 'LineWidth', 2.5, 'LineStyle', '--')
    colororder(line_color)
    axis([0 tSPAN(end) 0 inf]);
    title('Plot of cell populations - Adjusted Regimen', 'FontSize', 24);
    xlabel('time (hours)', 'FontSize', 18);
    ylabel('cell number', 'FontSize', 18);
    legend('susceptible', 'B resist', 'A resist', 'both resist', 'total');
    set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
    set(gca, 'FontSize', 16)
    hold off
    %save figure
    strFileName3 = sprintf('%s%s.fig', practicalSave, "_cellPopulationApproxTx");
    strFileName4 = sprintf('%s%s.png', practicalSave, "_cellPopulationApproxTx");
    saveas(gcf2,strFileName3); %saves population graph as fig
    saveas(gcf2,strFileName4); %saves population graph as png

end
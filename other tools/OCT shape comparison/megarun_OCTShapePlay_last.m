function megarun_OCTShapePlay_last()
    uMaxArr = [22.5, 45, 75, 150; 15, 30, 50, 100; 6, 12, 20, 40];
    [rowuMax, coluMax] = size(uMaxArr);
    threshVec = [0.25, 0.5, 0.8];

    comparisonRow = 7 * rowuMax;
    gammaEqComparisonArray = zeros(comparisonRow, coluMax);
    idealsEq = zeros(1, length(threshVec));
    worstsEq = zeros(1, length(threshVec));
    
    folderNameBase = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/manuscript figure generation/control shape comparison/last threshold break/"; %finish with a / so that the multiruns work
   
     for j=1:length(threshVec)
       for i=1:coluMax
            jointEqName = sprintf("joint equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqEqName = sprintf("sequential equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            figureName = sprintf("figures/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            mkdir(folderNameBase, jointEqName)
            mkdir(folderNameBase, seqEqName)
            mkdir(folderNameBase, figureName)
       end
     end


     for j=1:length(threshVec)
        for i=1:coluMax
            disp("this is run")
            disp((j-1)*(coluMax)+i)

            multisimRun_24_0828_jointWithSwitch_heatmapSave_GammaEqual_last(uMaxArr(j, i), threshVec(j), folderNameBase);
            multirun_24_0924_seqPopSwitch_gammaEqual_last(uMaxArr(j, i), threshVec(j), folderNameBase);
           
            %pull the population/time arrays for the best & worst of each run (or
            %maybe just best for now?)
        
            jointEqualHeatmapFile = sprintf("%sjoint equal/umax_%s_thresh_%s/simHeatmap.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqEqualHeatmapFile = sprintf("%ssequential equal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
        
            load(jointEqualHeatmapFile, 'heatmapVals');
            jointEqHeatmapVals = heatmapVals;
            load(seqEqualHeatmapFile, 'heatvectorVals');
            seqEqHeatmapVals = heatvectorVals;

            [jointEqMax, jointEqIndices] = max(jointEqHeatmapVals,[],"all");
            [rowJointEq, colJointEq] = ind2sub(size(jointEqHeatmapVals), jointEqIndices);

            close all

            %run the max version?
            [failureTime, switchpoint] = multirun_maximalOCTShape_last(uMaxArr(j, i), threshVec(j));

            %identify things to graph (rn 5 - joint best/worst, seq best/worst, max/max best)
                gammaEqComparisonArray(j, i) = jointEqMax;
                gammaEqComparisonArray((2*(length(threshVec)) + j), i) = min(nonzeros(jointEqHeatmapVals), [], "all");
                gammaEqComparisonArray((1*(length(threshVec)) + j), i) = max(seqEqHeatmapVals, [], "all");
                gammaEqComparisonArray((3*(length(threshVec)) + j), i) = min(nonzeros(seqEqHeatmapVals), [], "all");
                gammaEqComparisonArray((4*(length(threshVec)) + j), i) = multirun_uMaxA_uSingB_last(uMaxArr(j, i), threshVec(j));
                gammaEqComparisonArray((5*(length(threshVec)) + j), i) = multirun_elsaEqn(uMaxArr(j, i), threshVec(j));
                gammaEqComparisonArray((6*(length(threshVec)) + j), i) = failureTime;

            %calculate ideals/extremes
                %load population structs
                    jointEqualParamFile = sprintf("%sjoint equal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(jointEqualParamFile, 'ratio_2naught_3naught')
                    jointEqPopVec = ratio_2naught_3naught;
                    load(jointEqualParamFile, 'ratio_concA_concB')
                    jointEqConcVec = ratio_concA_concB;

                    jointEqConc = jointEqConcVec(rowJointEq);
                    jointEqPop = jointEqPopVec(colJointEq);

                    jointEqPopulationStruct = sprintf("%sjoint equal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_concAconcBratio%s_gammaEqual.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointEqPop), num2str(jointEqConc));
                    
                    load(jointEqPopulationStruct, 'popStruct');
                    jointEqP = popStruct;

                    idealEqual = calculateIdealFromPopStruct(jointEqP);

                    worstEqual = calculateWorstFromPopStruct(jointEqP);

                    idealsEq(j) = idealEqual;

                    worstsEq(j) = worstEqual;
        end
     end

    %graph the whole dang thing heheheheheheheh
close all
 
    figure(1)
    hold on;
        xAxisLabelNames = {};
        legendNames = {};
        for j=1:length(threshVec)
            for i=1:coluMax
                xLabelHere = {sprintf("threshold - %s, uMax - %s", num2str(threshVec(j)), num2str(uMaxArr(j, i)))};
                xAxisLabelNames = [xAxisLabelNames, xLabelHere];
                %graph the best and worst lines
                xaxisNum = i + (j - 1) * (2 + coluMax);
                xaxisLeft = xaxisNum - 0.07;
                xaxisRight = xaxisNum + 0.07;

                idealsHold = idealsEq(j);
                worstsHold = worstsEq(j);
                line([xaxisLeft, xaxisRight], [idealsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisLeft, xaxisRight], [worstsHold, worstsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisNum, xaxisNum], [worstsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off') 
            
                %graph the dots at that line
                plot((xaxisNum - .06), gammaEqComparisonArray(j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', 	"#0072BD")
                plot((xaxisNum - .04), gammaEqComparisonArray(1*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#4DBEEE")
                plot((xaxisNum - .02), gammaEqComparisonArray(2*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#3d630a")
                plot((xaxisNum), gammaEqComparisonArray(3*(length(threshVec))+j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#77AC30")
                plot((xaxisNum + .02), gammaEqComparisonArray(5*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#EDB120")
                plot((xaxisNum + .04), gammaEqComparisonArray(4*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#7E2F8E')
                plot((xaxisNum + .06), gammaEqComparisonArray(6*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#9a65b5')
            end
        end
        %axis settings
        maxIdeal = max(idealsEq, [], "all")+50;
        axis([0 18 0 maxIdeal])
        legendNames = [legendNames, {'Best Joint', 'Best Sequential', 'Worst Joint', 'Worst Sequential', 'Elsa Eqn', 'uMax/uSing->Max', 'Max/Max Opt'}];
        legend(legendNames)
        legend("Location", "northwest")

        legend show

        ylabel("time (hours)", 'FontSize', 18)
        xticks([1 2 3 4 7 8 9 10 13 14 15 16])
        xticklabels(xAxisLabelNames);

        title("Equal Gamma, All uMax & Thresh FIRST")

        set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
        set(gca, 'FontSize', 16)
        hold off
        equalGammaComparisonSave = sprintf("%sfigures/EqualGammaRegimenComparisonFigure", folderNameBase);
        %save figure
        strFileNameEqGCompFig = sprintf('%s.fig', equalGammaComparisonSave);
        strFileNameEqGCompPng = sprintf('%s.png', equalGammaComparisonSave);
        saveas(gcf,strFileNameEqGCompFig); %saves population graph as fig
        saveas(gcf,strFileNameEqGCompPng); %saves population graph as png
        %graph all of the above & save?

end
function megarun_OCTShapePlay_AandB_EqandUneq()
    uMaxArr = [22.5, 45, 75, 150; 15, 30, 50, 100; 6, 12, 20, 40];
    %uMaxArr = [6, 12, 20, 40];
    %uMaxArr = [40];
    [rowuMax, coluMax] = size(uMaxArr);
    threshVec = [0.25, 0.5, 0.8];
    %threshVec = [0.8];

    comparisonRow = 11 * rowuMax;
    gammaEqComparisonArray = zeros(comparisonRow, coluMax);
    gammaUneqComparisonArray = zeros(comparisonRow, coluMax);
    idealsEq = zeros(1, length(threshVec));
    worstsEq = zeros(1, length(threshVec));
    idealsUneq = zeros(1, length(threshVec));
    worstsUneq = zeros(1, length(threshVec));
    
    folderNameBase = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/manuscript/oct regimen comparison/1015 full run/"; %finish with a / so that the multiruns work
   
     for j=1:length(threshVec)
       for i=1:coluMax
            jointEqName = sprintf("joint equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            jointUneqName = sprintf("joint unequal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqEqName = sprintf("sequential equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqUneqName = sprintf("sequential unequal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqBEqName = sprintf("sequential b equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqBUneqName = sprintf("sequential b unequal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            figureName = sprintf("figures/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            maxmaxAName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualMaxMaxA", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            maxmaxBName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualMaxMaxB", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            maxmaxAUneqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualMaxMaxA_uneq", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            maxmaxBUneqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualMaxMaxB_uneq", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            elsaEqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualElsaEqual", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            elsaUneqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualElsaUnequal", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            uMaxAuSingBName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualuMaxAuSingB", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            uMaxAuSingBUneqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualuMaxAuSingBuneq", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            uMaxBuSingAName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualuMaxBuSingA", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            uMaxBuSingAUneqName = sprintf("otherOCTshapes/umax_%s_thresh_%s/individualuMaxBuSingAuneq", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            
            mkdir(folderNameBase, jointEqName)
            mkdir(folderNameBase, jointUneqName)
            mkdir(folderNameBase, seqEqName)
            mkdir(folderNameBase, seqUneqName)
            mkdir(folderNameBase, seqBEqName)
            mkdir(folderNameBase, seqBUneqName)
            mkdir(folderNameBase, figureName)
            mkdir(folderNameBase, maxmaxAName)
            mkdir(folderNameBase, maxmaxBName)
            mkdir(folderNameBase, maxmaxAUneqName)
            mkdir(folderNameBase, maxmaxBUneqName)
            mkdir(folderNameBase, elsaEqName)
            mkdir(folderNameBase, elsaUneqName)
            mkdir(folderNameBase, uMaxAuSingBName)
            mkdir(folderNameBase, uMaxAuSingBUneqName)
            mkdir(folderNameBase, uMaxBuSingAName)
            mkdir(folderNameBase, uMaxBuSingAUneqName)
       end
     end


     for j=1:length(threshVec)
        for i=1:coluMax
            disp("this is run")
            disp((j-1)*(coluMax)+i)

            multisimRun_24_0617_jointWithSwitch_heatmapSave_GammaEqual(uMaxArr(j, i), threshVec(j), folderNameBase);
            multisimRun_24_0617_jointWithSwitch_heatmapSave_GammaUnequal(uMaxArr(j, i), threshVec(j), folderNameBase);
         
            multirun_24_0918_seqPopSwitch_gammaEqual(uMaxArr(j, i), threshVec(j), folderNameBase);
            multirun_24_0918_seqPopSwitch_gammaUnequal(uMaxArr(j, i), threshVec(j), folderNameBase);
            %disp("normal Seq Pop")
            multirun_24_1002_seqPopSwitchB_gammaEqual(uMaxArr(j, i), threshVec(j), folderNameBase);
            multirun_24_0918_seqPopSwitchB_gammaUnequal(uMaxArr(j, i), threshVec(j), folderNameBase);
            %disp("normal Seq B Pop")

            %pull the population/time arrays for the best & worst of each run (or
            %maybe just best for now?)
        
           jointEqualHeatmapFile = sprintf("%sjoint equal/umax_%s_thresh_%s/simHeatmap.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           jointUnequalHeatmapFile = sprintf("%sjoint unequal/umax_%s_thresh_%s/simHeatmap.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqEqualHeatmapFile = sprintf("%ssequential equal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqUnequalHeatmapFile = sprintf("%ssequential unequal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqBEqualHeatmapFile = sprintf("%ssequential b equal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqBUnequalHeatmapFile = sprintf("%ssequential b unequal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
        

           load(jointEqualHeatmapFile, 'heatmapVals');
           jointEqHeatmapVals = heatmapVals;
           load(jointUnequalHeatmapFile, 'heatmapVals');
           jointUneqHeatmapVals = heatmapVals;
           load(seqEqualHeatmapFile, 'heatvectorVals');
           seqEqHeatmapVals = heatvectorVals;
           load(seqUnequalHeatmapFile, 'heatvectorVals');
           seqUneqHeatmapVals = heatvectorVals;
           load(seqBEqualHeatmapFile, 'heatvectorVals');
           seqBEqHeatmapVals = heatvectorVals;
           load(seqBUnequalHeatmapFile, 'heatvectorVals');
           seqBUneqHeatmapVals = heatvectorVals;

           %find best option
           [jointEqMax, jointEqIndices] = max(jointEqHeatmapVals,[],"all");
           [rowJointEq, colJointEq] = ind2sub(size(jointEqHeatmapVals), jointEqIndices);
           [jointUneqMax, jointUneqIndices] = max(jointUneqHeatmapVals,[],"all");
           [rowJointUneq, colJointUneq] = ind2sub(size(jointUneqHeatmapVals), jointUneqIndices);
           [seqEqMax, seqEqIndices] = max(seqEqHeatmapVals,[],"all");
           [rowSeqEq, colSeqEq] = ind2sub(size(seqEqHeatmapVals), seqEqIndices);
           [seqUneqMax, seqUneqIndices] = max(seqUneqHeatmapVals,[],"all");
           [rowSeqUneq, colSeqUneq] = ind2sub(size(seqUneqHeatmapVals), seqUneqIndices);
           [seqBEqMax, seqBEqIndices] = max(seqBEqHeatmapVals,[],"all");
           [rowSeqBEq, colSeqBEq] = ind2sub(size(seqBEqHeatmapVals), seqBEqIndices);
           [seqBUneqMax, seqBUneqIndices] = max(seqBUneqHeatmapVals,[],"all");
           [rowSeqBUneq, colSeqBUneq] = ind2sub(size(seqBUneqHeatmapVals), seqBUneqIndices);

            
            close all

            %run the max version?
            [failureTimeA, switchpoint] = multirun_maximalOCTShape(uMaxArr(j, i), threshVec(j), folderNameBase);
            [failureTimeB, switchpoint2] = multirun_maximalOCTShape_Bfirst(uMaxArr(j, i), threshVec(j), folderNameBase);
            [uneqFailureTimeA, switchpointu] = multirun_maximalOCTShape_uneq(uMaxArr(j, i), threshVec(j), folderNameBase);
            [uneqFailureTimeB, switchpoint2u] = multirun_maximalOCTShape_Bfirst_uneq(uMaxArr(j, i), threshVec(j), folderNameBase);
            

            %identify things to graph (rn 5 - joint best/worst, seq best/worst, max/max best)
                gammaEqComparisonArray(j, i) = jointEqMax;
                gammaEqComparisonArray((1*(length(threshVec)) + j), i) = seqEqMax;
                gammaEqComparisonArray((2*(length(threshVec)) + j), i) = seqBEqMax;
                gammaEqComparisonArray((3*(length(threshVec)) + j), i) = min(nonzeros(jointEqHeatmapVals), [], "all");
                gammaEqComparisonArray((4*(length(threshVec)) + j), i) = min(nonzeros(seqEqHeatmapVals), [], "all");
                gammaEqComparisonArray((5*(length(threshVec)) + j), i) = min(nonzeros(seqBEqHeatmapVals), [], "all");
                gammaEqComparisonArray((6*(length(threshVec)) + j), i) = multirun_elsaEqn(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaEqComparisonArray((7*(length(threshVec)) + j), i) = multirun_uMaxA_uSingB(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaEqComparisonArray((8*(length(threshVec)) + j), i) = multirun_uMaxB_uSingA(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaEqComparisonArray((9*(length(threshVec)) + j), i) = failureTimeA;
                gammaEqComparisonArray((10*(length(threshVec)) + j), i) = failureTimeB;

                gammaUneqComparisonArray(j, i) = jointUneqMax;
                gammaUneqComparisonArray((1*(length(threshVec)) + j), i) = seqUneqMax;
                gammaUneqComparisonArray((2*(length(threshVec)) + j), i) = seqBUneqMax;
                gammaUneqComparisonArray((3*(length(threshVec)) + j), i) = min(nonzeros(jointUneqHeatmapVals), [], "all");
                gammaUneqComparisonArray((4*(length(threshVec)) + j), i) = min(nonzeros(seqUneqHeatmapVals), [], "all");
                gammaUneqComparisonArray((5*(length(threshVec)) + j), i) = min(nonzeros(seqBUneqHeatmapVals), [], "all");
                gammaUneqComparisonArray((6*(length(threshVec)) + j), i) = multirun_elsaEqn_uneq(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaUneqComparisonArray((7*(length(threshVec)) + j), i) = multirun_uMaxA_uSingB_uneq(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaUneqComparisonArray((8*(length(threshVec)) + j), i) = multirun_uMaxB_uSingA_uneq(uMaxArr(j, i), threshVec(j), folderNameBase);
                gammaUneqComparisonArray((9*(length(threshVec)) + j), i) = uneqFailureTimeA;
                gammaUneqComparisonArray((10*(length(threshVec)) + j), i) = uneqFailureTimeB;

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

                    jointUnequalParamFile = sprintf("%sjoint unequal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(jointUnequalParamFile, 'ratio_2naught_3naught')
                    jointUneqPopVec = ratio_2naught_3naught;
                    load(jointUnequalParamFile, 'ratio_concA_concB')
                    jointUneqConcVec = ratio_concA_concB;

                    jointUneqConc = jointUneqConcVec(rowJointUneq);
                    jointUneqPop = jointUneqPopVec(colJointUneq);

                    jointUneqPopulationStruct = sprintf("%sjoint unequal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_concAconcBratio%s_gammaUnequal.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointEqPop), num2str(jointEqConc));
                    
                    load(jointUneqPopulationStruct, 'popStruct');
                    jointUneqP = popStruct;

                    idealUnequal = calculateIdealFromPopStruct(jointUneqP);

                    worstUnequal = calculateWorstFromPopStruct(jointUneqP);

                    idealsUneq(j) = idealUnequal;

                    worstsUneq(j) = worstUnequal;
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
                xaxisLeft = xaxisNum - 0.1;
                xaxisRight = xaxisNum + 0.1;

                idealsHold = idealsEq(j);
                worstsHold = worstsEq(j);
                line([xaxisLeft, xaxisRight], [idealsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisLeft, xaxisRight], [worstsHold, worstsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisNum, xaxisNum], [worstsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off') 
            
                %graph the dots at that line
                plot((xaxisNum - .1), gammaEqComparisonArray(j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', 	"#0072BD")
                plot((xaxisNum - .08), gammaEqComparisonArray(1*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#4DBEEE")
                plot((xaxisNum - .06), gammaEqComparisonArray(2*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', 	"#4DBEEE", 'MarkerFaceColor', "#4DBEEE")
                plot((xaxisNum - .04), gammaEqComparisonArray(3*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#3d630a")
                plot((xaxisNum - .02), gammaEqComparisonArray(4*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#77AC30")
                plot((xaxisNum), gammaEqComparisonArray(5*(length(threshVec))+j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', "#77AC30", 'MarkerFaceColor', "#77AC30")
                plot((xaxisNum + .02), gammaEqComparisonArray(6*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#EDB120")
                plot((xaxisNum + .04), gammaEqComparisonArray(7*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#9a65b5')
                plot((xaxisNum + .06), gammaEqComparisonArray(8*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', '#9a65b5', 'MarkerFaceColor', "#9a65b5")
                plot((xaxisNum + .08), gammaEqComparisonArray(9*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#7E2F8E')
                plot((xaxisNum + .1), gammaEqComparisonArray(10*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', '#7E2F8E', 'MarkerFaceColor', "#7E2F8E")
            end
        end
        %axis settings
        maxIdeal = max(idealsEq, [], "all")+50;
        axis([0 18 0 maxIdeal])
        legendNames = [legendNames, {'Best Joint', 'Best Sequential', 'Best Sequential B', 'Worst Joint', 'Worst Sequential', 'Worst Sequential B', 'Elsa Eqn', 'uMax/uSing->Max', 'uSing->Max/uMax', 'Max/Max Opt A first', 'Max/Max Opt B first'}];
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
        savefig(gcf,strFileNameEqGCompFig); %saves population graph as fig
        %saveas(gcf,strFileNameEqGCompPng); %saves population graph as png
        %graph all of the above & save?


    figure(2)
    hold on;
        xAxisLabelNames = {};
        legendNames = {};
        for j=1:length(threshVec)
            for i=1:coluMax
                xLabelHere = {sprintf("threshold - %s, uMax - %s", num2str(threshVec(j)), num2str(uMaxArr(j, i)))};
                xAxisLabelNames = [xAxisLabelNames, xLabelHere];
                %graph the best and worst lines
                xaxisNum = i + (j - 1) * (2 + coluMax);
                xaxisLeft = xaxisNum - 0.1;
                xaxisRight = xaxisNum + 0.1;

                idealsHold = idealsUneq(j);
                worstsHold = worstsUneq(j);
                line([xaxisLeft, xaxisRight], [idealsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisLeft, xaxisRight], [worstsHold, worstsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisNum, xaxisNum], [worstsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off') 
            
                %graph the dots at that line
                plot((xaxisNum - .1), gammaUneqComparisonArray(j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', 	"#0072BD")
                plot((xaxisNum - .08), gammaUneqComparisonArray(1*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#4DBEEE")
                plot((xaxisNum - .06), gammaUneqComparisonArray(2*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', 	"#4DBEEE", 'MarkerFaceColor', "#4DBEEE")
                plot((xaxisNum - .04), gammaUneqComparisonArray(3*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#3d630a")
                plot((xaxisNum - .02), gammaUneqComparisonArray(4*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#77AC30")
                plot((xaxisNum), gammaUneqComparisonArray(5*(length(threshVec))+j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', "#77AC30", 'MarkerFaceColor', "#77AC30")
                plot((xaxisNum + .02), gammaUneqComparisonArray(6*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#EDB120")
                plot((xaxisNum + .04), gammaUneqComparisonArray(7*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#9a65b5')
                plot((xaxisNum + .06), gammaUneqComparisonArray(8*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', '#9a65b5', 'MarkerFaceColor', "#9a65b5")
                plot((xaxisNum + .08), gammaUneqComparisonArray(9*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#7E2F8E')
                plot((xaxisNum + .1), gammaUneqComparisonArray(10*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '^', 'MarkerSize', 10, 'Color', '#7E2F8E', 'MarkerFaceColor', "#7E2F8E")
            end
        end
        %axis settings
        maxIdeal = max(idealsUneq, [], "all")+50;
        axis([0 18 0 maxIdeal])
        legendNames = [legendNames, {'Best Joint', 'Best Sequential', 'Best Sequential B', 'Worst Joint', 'Worst Sequential', 'Worst Sequential B', 'Elsa Eqn', 'uMax/uSing->Max', 'uSing->Max/uMax', 'Max/Max Opt A first', 'Max/Max Opt B first'}];
        legend(legendNames)
        legend("Location", "northwest")

        legend show

        ylabel("time (hours)", 'FontSize', 18)
        xticks([1 2 3 4 7 8 9 10 13 14 15 16])
        xticklabels(xAxisLabelNames);

        title("Unequal Gamma, All uMax & Thresh FIRST")

        set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
        set(gca, 'FontSize', 16)
        hold off
        unequalGammaComparisonSave = sprintf("%sfigures/UnequalGammaRegimenComparisonFigure", folderNameBase);
        %save figure
        strFileNameUneqGCompFig = sprintf('%s.fig', unequalGammaComparisonSave);
        strFileNameUneqGCompPng = sprintf('%s.png', unequalGammaComparisonSave);
        savefig(gcf,strFileNameUneqGCompFig); %saves population graph as fig
        %saveas(gcf,strFileNameUneqGCompPng); %saves population graph as png
        %graph all of the above & save?

end
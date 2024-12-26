function megarun_practicality_24_0816_uMaxkR_fixed()
    %uMaxVec = [20, 40, 80];
    uMaxArr = [22.5, 45, 75, 150; 15, 30, 50, 100; 6, 12, 20, 40];
    [rowuMax, coluMax] = size(uMaxArr);
    threshVec = [0.25, 0.5, 0.8];
    timeStep = 24; %in hours
    dosingIncrement = 2.5; %in nM

%    uMaxVec = [0.5, 0.75, 1];
%    threshVec = [0.25, 0.5, 0.8];

    comparisonRow = 8 * rowuMax;
    gammaEqComparisonArray = zeros(comparisonRow, coluMax);
    gammaUneqComparisonArray = zeros(comparisonRow, coluMax);
    idealsEq = zeros(1, length(threshVec));
    idealsUneq = zeros(1, length(threshVec));
    worstsEq = zeros(1, length(threshVec));
    worstsUneq = zeros(1, length(threshVec));

   folderNameBase = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/imlosinit/hey 3/"; %finish with a / so that the multiruns work
   
   for j=1:length(threshVec)
       for i=1:coluMax
            jointEqName = sprintf("joint equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            jointUneqName = sprintf("joint unequal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqEqName = sprintf("sequential equal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            seqUneqName = sprintf("sequential unequal/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            figureName = sprintf("figures/umax_%s_thresh_%s/individualGraphs", num2str(uMaxArr(j, i)), num2str(threshVec(j)));
            mkdir(folderNameBase, jointEqName)
            mkdir(folderNameBase, jointUneqName)
            mkdir(folderNameBase, seqEqName)
            mkdir(folderNameBase, seqUneqName)
            mkdir(folderNameBase, figureName)
       end
   end

   for j=1:length(threshVec)
       for i=1:coluMax
           disp("this is run")
           disp((j-1)*(coluMax)+i)

           multisimRun_24_0617_jointWithSwitch_heatmapSave_GammaEqual(uMaxArr(j, i), threshVec(j), folderNameBase);
           multisimRun_24_0617_jointWithSwitch_heatmapSave_GammaUnequal(uMaxArr(j, i), threshVec(j), folderNameBase);
           multirun_24_0617_seqBWithEpsilonNoSwitch_GammaEqual(uMaxArr(j, i), threshVec(j), folderNameBase);
           multirun_24_0617_seqBWithEpsilonNoSwitch_GammaUnequal(uMaxArr(j, i), threshVec(j), folderNameBase);
           
           %pull the population/time arrays for the best & worst of each run (or
           %maybe just best for now?)
        
           jointEqualHeatmapFile = sprintf("%sjoint equal/umax_%s_thresh_%s/simHeatmap.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           jointUnequalHeatmapFile = sprintf("%sjoint unequal/umax_%s_thresh_%s/simHeatmap.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqEqualHeatmapFile = sprintf("%ssequential equal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
           seqUnequalHeatmapFile = sprintf("%ssequential unequal/umax_%s_thresh_%s/simHeatvector.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
        
           load(jointEqualHeatmapFile, 'heatmapVals');
           jointEqHeatmapVals = heatmapVals;
           load(jointUnequalHeatmapFile, 'heatmapVals');
           jointUneqHeatmapVals = heatmapVals;
           load(seqEqualHeatmapFile, 'heatvectorVals');
           seqEqHeatmapVals = heatvectorVals;
           load(seqUnequalHeatmapFile, 'heatvectorVals');
           seqUneqHeatmapVals = heatvectorVals;

           %find best option
           [jointEqMax, jointEqIndices] = max(jointEqHeatmapVals,[],"all");
           [rowJointEq, colJointEq] = ind2sub(size(jointEqHeatmapVals), jointEqIndices);
           [jointUneqMax, jointUneqIndices] = max(jointUneqHeatmapVals,[],"all");
           [rowJointUneq, colJointUneq] = ind2sub(size(jointUneqHeatmapVals), jointUneqIndices);
           [seqEqMax, seqEqIndices] = max(seqEqHeatmapVals,[],"all");
           [rowSeqEq, colSeqEq] = ind2sub(size(seqEqHeatmapVals), seqEqIndices);
           [seqUneqMax, seqUneqIndices] = max(seqUneqHeatmapVals,[],"all");
           [rowSeqUneq, colSeqUneq] = ind2sub(size(seqUneqHeatmapVals), seqUneqIndices);

           %run the practical version of the best
                jointEqualParamFile = sprintf("%sjoint equal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(jointEqualParamFile, 'ratio_2naught_3naught')
                    jointEqPopVec = ratio_2naught_3naught;
                    load(jointEqualParamFile, 'ratio_concA_concB')
                    jointEqConcVec = ratio_concA_concB;
                jointUnequalParamFile = sprintf("%sjoint unequal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(jointUnequalParamFile, 'ratio_2naught_3naught')
                    jointUneqPopVec = ratio_2naught_3naught;
                    load(jointUnequalParamFile, 'ratio_concA_concB')
                    jointUneqConcVec = ratio_concA_concB;
                seqEqualParamFile = sprintf("%ssequential equal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(seqEqualParamFile, 'ratio_2naught_3naught')
                    seqEqPopVec = ratio_2naught_3naught;
                seqUnequalParamFile = sprintf("%ssequential unequal/umax_%s_thresh_%s/VectorRange.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    load(seqUnequalParamFile, 'ratio_2naught_3naught')
                    seqUneqPopVec = ratio_2naught_3naught;
                %open the file with the vectors of parameters & find the
                %mix that resulted in the best ** note, this is a
                %symmetrical thing now
                    %joint equal
                    jointEqConc = jointEqConcVec(rowJointEq);
                    jointEqPop = jointEqPopVec(colJointEq);

                    jointEqResultsFile = sprintf("%sjoint equal/umax_%s_thresh_%s/individualGraphs/test_x2x3ratio%s_concAconcBratio%s_gammaEqual_timeAndPopArr.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointEqPop), num2str(jointEqConc));
                    jointEqPopulationStruct = sprintf("%sjoint equal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_concAconcBratio%s_gammaEqual.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointEqPop), num2str(jointEqConc));
                    
                    jointEqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/jointEqualBestPractical_timestep_%sh_dosingInc_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(dosingIncrement));

                    %joint unequal
                    jointUneqConc = jointUneqConcVec(rowJointUneq);
                    jointUneqPop = jointUneqPopVec(colJointUneq);

                    jointUneqResultsFile = sprintf("%sjoint unequal/umax_%s_thresh_%s/individualGraphs/test_x2x3ratio%s_concAconcBratio%s_gammaUnequal_timeAndPopArr.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointUneqPop), num2str(jointUneqConc));
                    jointUneqPopulationStruct = sprintf("%sjoint unequal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_concAconcBratio%s_gammaUnequal.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(jointUneqPop), num2str(jointUneqConc));
                    
                    jointUneqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/jointUnequalBestPractical_timestep_%sh_dosingInc_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(dosingIncrement));

                    %sequential equal
                    seqEqPop = seqEqPopVec(rowSeqEq);

                    seqEqResultsFile = sprintf("%ssequential equal/umax_%s_thresh_%s/individualGraphs/test_x2x3ratio%s_gammaEqual_timeAndPopArr.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(seqEqPop));
                    seqEqPopulationStruct = sprintf("%ssequential equal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_gammaEqual.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(seqEqPop));
                   
                    seqEqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/seqEqualBestPractical_timestep_%sh_dosingInc_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(dosingIncrement));

                    %sequential unequal
                    seqUneqPop = seqUneqPopVec(rowSeqUneq);

                    seqUneqResultsFile = sprintf("%ssequential unequal/umax_%s_thresh_%s/individualGraphs/test_x2x3ratio%s_gammaUnequal_timeAndPopArr.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(seqUneqPop));
                    seqUneqPopulationStruct = sprintf("%ssequential unequal/umax_%s_thresh_%s/individualGraphs/population_x2x3ratio%s_gammaUnequal.mat", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(seqUneqPop));
               
                    seqUneqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/seqUnequalBestPractical_timestep_%sh_dosingInc_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(dosingIncrement));

                %feed population & time array into the joint or the
                %sequential practicality function

                    jointEqThresh = practicalityAnalysisJoint(jointEqPopulationStruct, jointEqResultsFile, jointEqSaveFile, timeStep, jointEqConc, dosingIncrement);
                    jointUneqThresh = practicalityAnalysisJoint(jointUneqPopulationStruct, jointUneqResultsFile, jointUneqSaveFile, timeStep, jointUneqConc, dosingIncrement);
                    seqEqThresh = practicalityAnalysisSeq(seqEqPopulationStruct, seqEqResultsFile, seqEqSaveFile, timeStep, dosingIncrement);
                    seqUneqThresh = practicalityAnalysisSeq(seqUneqPopulationStruct, seqUneqResultsFile, seqUneqSaveFile, timeStep, dosingIncrement);

               %generate maximal treatment version
                    maxDoseA = 400;
                    maxDoseB = 400;
                    eqMaxComboSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/equalMaxCombo_timestep_%sh_doseA_%snM_doseB_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(maxDoseA), num2str(maxDoseB));
                    uneqMaxComboSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/unequalMaxCombo_timestep_%sh_doseA_%snM_doseB_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(maxDoseA), num2str(maxDoseB));
                    eqMaxSeqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/equalMaxSeq_timestep_%sh_doseA_%snM_doseB_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(maxDoseA), num2str(maxDoseB));
                    uneqMaxSeqSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/unequalMaxSeq_timestep_%sh_doseA_%snM_doseB_%snM", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)), num2str(timeStep), num2str(maxDoseA), num2str(maxDoseB));

                    
                    eqMaxComboThresh = maximalTreatmentCombo(jointEqPopulationStruct, jointEqResultsFile, eqMaxComboSaveFile, maxDoseA, maxDoseB, timeStep);
                    uneqMaxComboThresh = maximalTreatmentCombo(jointUneqPopulationStruct, jointUneqResultsFile, uneqMaxComboSaveFile, maxDoseA, maxDoseB, timeStep);
                    eqMaxSeqThresh = maximalTreatmentSequential(jointEqPopulationStruct, jointEqResultsFile, eqMaxSeqSaveFile, maxDoseA, maxDoseB, timeStep);
                    uneqMaxSeqThresh = maximalTreatmentSequential(jointUneqPopulationStruct, jointUneqResultsFile, uneqMaxSeqSaveFile, maxDoseA, maxDoseB, timeStep);

                %save figures & timepoint generated so that it can be fed
                %into final figure generation down below vv
                    %identify extremes
                        %load population structs
                            load(jointEqPopulationStruct, 'popStruct');
                            jointEqP = popStruct;
                            load(jointUneqPopulationStruct, 'popStruct');
                            jointUneqP = popStruct;

                            idealEqual = calculateIdealFromPopStruct(jointEqP);
                            idealUnequal = calculateIdealFromPopStruct(jointUneqP);

                            worstEqual = calculateWorstFromPopStruct(jointEqP);
                            worstUnequal = calculateWorstFromPopStruct(jointUneqP);

                            idealsEq(j) = idealEqual;
                            idealsUneq(j) = idealUnequal;

                            worstsEq(j) = worstEqual;
                            worstsUneq(j) = worstUnequal;

                    %collates the heatmap values
                    jointVec = [max(jointEqHeatmapVals, [], "all"), max(jointUneqHeatmapVals, [], "all")];
                    jointPracticalVec = [jointEqThresh, jointUneqThresh];
                    jointWorstVec = [min(nonzeros(jointEqHeatmapVals), [], "all"), min(nonzeros(jointUneqHeatmapVals), [], "all")];
                    seqVec = [max(seqEqHeatmapVals, [], "all"), max(seqUneqHeatmapVals, [], "all")];
                    seqPracticalVec = [seqEqThresh, seqUneqThresh];
                    seqWorstVec = [min(nonzeros(seqEqHeatmapVals), [], "all"), min(nonzeros(seqUneqHeatmapVals), [], "all")];
                    maximalComboVec = [eqMaxComboThresh, uneqMaxComboThresh];
                    maximalSeqVec = [eqMaxSeqThresh, uneqMaxSeqThresh];

                    close all

                    %make vectors
                    gammaEqComparisonArray(j, i) = jointVec(1);
                    gammaUneqComparisonArray(j, i) = jointVec(2);

                    gammaEqComparisonArray((1*(length(threshVec)) + j), i) = seqVec(1);
                    gammaUneqComparisonArray((1*(length(threshVec)) + j), i) = seqVec(2);

                    gammaEqComparisonArray((2*(length(threshVec)) + j), i) = jointWorstVec(1);
                    gammaUneqComparisonArray((2*(length(threshVec)) + j), i) = jointWorstVec(2);

                    gammaEqComparisonArray((3*(length(threshVec)) + j), i) = seqWorstVec(1);
                    gammaUneqComparisonArray((3*(length(threshVec)) + j), i) = seqWorstVec(2);

                    gammaEqComparisonArray((4*(length(threshVec)) + j), i) = jointEqThresh;
                    gammaUneqComparisonArray((4*(length(threshVec)) + j), i) = jointUneqThresh;

                    gammaEqComparisonArray((5*(length(threshVec)) + j), i) = seqEqThresh;
                    gammaUneqComparisonArray((5*(length(threshVec)) + j), i) = seqUneqThresh;

                    gammaEqComparisonArray((6*(length(threshVec)) + j), i) = eqMaxComboThresh;
                    gammaUneqComparisonArray((6*(length(threshVec)) + j), i) = uneqMaxComboThresh;

                    gammaEqComparisonArray((7*(length(threshVec)) + j), i) = eqMaxSeqThresh;
                    gammaUneqComparisonArray((7*(length(threshVec)) + j), i) = uneqMaxSeqThresh;


                    comparisonSaveFile = sprintf("%sfigures/umax_%s_thresh_%s/RegimenComparisonFigure", folderNameBase, num2str(uMaxArr(j, i)), num2str(threshVec(j)));
                    
                
               figure((j-1)*length(threshVec)+i)
                hold on

                %best & worst lines
                line([0.925, 1.075], [idealEqual, idealEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([0.925, 1.075], [worstEqual, worstEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1, 1], [worstEqual, idealEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1.925, 2.075], [idealUnequal, idealUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1.925, 2.075], [worstUnequal, worstUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([2, 2], [worstUnequal, idealUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')

                %plot the various bits and bobs
                plot([0.925, 1.925], jointVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Best Joint')
                plot([0.946, 1.946], seqVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Best Sequential')
                plot([0.968, 1.968], jointWorstVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Worst Joint')
                plot([0.989, 1.989], seqWorstVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Worst Sequential')
                plot([1.011, 2.011], jointPracticalVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Practical Joint')
                plot([1.032, 2.032], seqPracticalVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Practical Sequential')
                plot([1.054, 2.054], maximalSeqVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Maximal Sequential')
                plot([1.075, 2.075], maximalComboVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Maximal Combination', 'Color', '#262351')

                %axis settings
                maxIdeal = max(idealUnequal, idealEqual)+50;
                axis([0 3 0 maxIdeal])

                legend show

                ylabel("time (hours)", 'FontSize', 18)
                xticks([1 2])
                xticklabels({"Equal \gamma", "Unequal \gamma"});

                title(sprintf("uMax - %s Threshold - %s",num2str(uMaxArr(j, i)), num2str(threshVec(j))))

                set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
                set(gca, 'FontSize', 16)
                hold off
                %save figure
                strFileNameCompFig = sprintf('%s.fig', comparisonSaveFile);
                strFileNameCompPng = sprintf('%s.png', comparisonSaveFile);
                saveas(gcf,strFileNameCompFig); %saves population graph as fig
                saveas(gcf,strFileNameCompPng); %saves population graph as png

       end
   end

   figure(15)
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
            
%                 line([xaxisLeft, xaxisRight], [idealsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([xaxisLeft, xaxisRight], [worstsEq(j), worstsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 
%                 line([(xaxisNum - 0.07), (xaxisNum + 0.07)], [idealsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([(xaxisNum - 0.07), (xaxisNum + 0.07)], [worstsEq(j), worstsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([xaxisNum, xaxisNum], [worstsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')   
            
                %graph the dots at that line
                plot((xaxisNum - .07), gammaEqComparisonArray(j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', 	"#0072BD")
                plot((xaxisNum - .05), gammaEqComparisonArray(length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#4DBEEE")
                plot((xaxisNum - .03), gammaEqComparisonArray(2*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#3d630a")
                plot((xaxisNum - .01), gammaEqComparisonArray(3*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#77AC30")
                plot((xaxisNum + .01), gammaEqComparisonArray(4*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#d8861a")
                plot((xaxisNum + .03), gammaEqComparisonArray(5*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#EDB120')
                plot((xaxisNum + .05), gammaEqComparisonArray(6*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#7E2F8E')
                plot((xaxisNum + .07), gammaEqComparisonArray(7*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#9a65b5')
            end
        end
        %axis settings
        maxIdeal = max(idealsEq, [], "all")+50;
        axis([0 18 0 maxIdeal])
        legendNames = [legendNames, {'Best Joint', 'Best Sequential', 'Worst Joint', 'Worst Sequential', 'Practical Joint', 'Practical Sequential', 'Maximal Combo', 'Maximal Sequential'}];
        legend(legendNames)
        legend("Location", "northwest")

        legend show

        ylabel("time (hours)", 'FontSize', 18)
        xticks([1 2 3 4 7 8 9 10 13 14 15 16])
        xticklabels(xAxisLabelNames);

        title("Equal Gamma, All uMax & Thresh")

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


   figure(16)
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

                idealsHold = idealsUneq(j);
                worstsHold = worstsUneq(j);
                line([xaxisLeft, xaxisRight], [idealsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisLeft, xaxisRight], [worstsHold, worstsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([xaxisNum, xaxisNum], [worstsHold, idealsHold], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')   
            
%                 line([xaxisLeft, xaxisRight], [idealsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([xaxisLeft, xaxisRight], [worstsEq(j), worstsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 
%                 line([(xaxisNum - 0.07), (xaxisNum + 0.07)], [idealsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([(xaxisNum - 0.07), (xaxisNum + 0.07)], [worstsEq(j), worstsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
%                 line([xaxisNum, xaxisNum], [worstsEq(j), idealsEq(j)], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')   
            
                %graph the dots at that line
                plot((xaxisNum - .07), gammaUneqComparisonArray(j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', 	"#0072BD")
                plot((xaxisNum - .05), gammaUneqComparisonArray(length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#4DBEEE")
                plot((xaxisNum - .03), gammaUneqComparisonArray(2*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#3d630a")
                plot((xaxisNum - .01), gammaUneqComparisonArray(3*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#77AC30")
                plot((xaxisNum + .01), gammaUneqComparisonArray(4*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', "#d8861a")
                plot((xaxisNum + .03), gammaUneqComparisonArray(5*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#EDB120')
                plot((xaxisNum + .05), gammaUneqComparisonArray(6*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#7E2F8E')
                plot((xaxisNum + .07), gammaUneqComparisonArray(7*length(threshVec) + j, i), 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'Color', '#9a65b5')
            end
        end
        %axis settings
        maxIdeal = max(idealsUneq, [], "all")+50;
        axis([0 18 0 maxIdeal])
        legendNames = [legendNames, {'Best Joint', 'Best Sequential', 'Worst Joint', 'Worst Sequential', 'Practical Joint', 'Practical Sequential', 'Maximal Combo', 'Maximal Sequential'}];
        legend(legendNames)
        legend("Location", "northwest")

        legend show

        ylabel("time (hours)", 'FontSize', 18)
        xticks([1 2 3 4 7 8 9 10 13 14 15 16])
        xticklabels(xAxisLabelNames);

        title("Unequal Gamma, All uMax & Thresh")

        set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
        set(gca, 'FontSize', 16)
        hold off
        unequalGammaComparisonSave = sprintf("%sfigures/UnequalGammaRegimenComparisonFigure", folderNameBase);
        %save figure
        strFileNameUneqGCompFig = sprintf('%s.fig', unequalGammaComparisonSave);
        strFileNameUneqGCompPng = sprintf('%s.png', unequalGammaComparisonSave);
        saveas(gcf,strFileNameUneqGCompFig); %saves population graph as fig
        saveas(gcf,strFileNameUneqGCompPng); %saves population graph as png
        %graph all of the above & save?
end
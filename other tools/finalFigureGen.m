function finalFigureGen()

    %base filenames
    folderJointEqual = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/final round threshold error rerun/joint equal/";
    folderJointUnequal = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/final round threshold error rerun/joint unequal/";
    folderSeqEqual = "/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/final round threshold error rerun/sequential equal/";
    folderSeqUnequal ="/Users/aftonwiddershins/Desktop/stuff/stuff/academic/phd/thesis/Mirror/thesis data/model data/double bound generation/final round threshold error rerun/sequential unequal/";

    %umax & thresh variables
    uMaxVec = [0.5, 0.75, 1];
    threshVec = [0.25, 0.5, 0.8];

    
    
    %base figure
    for i=1:length(uMaxVec)
        for j=1:length(threshVec)
            %full file names
            fileJointEqual = sprintf("%sumax_%s_thresh_%s/", folderJointEqual, num2str(uMaxVec(i)), num2str(threshVec(j)));
            fileJointUnequal = sprintf("%sumax_%s_thresh_%s/", folderJointUnequal, num2str(uMaxVec(i)), num2str(threshVec(j)));
            fileSeqEqual = sprintf("%sumax_%s_thresh_%s/", folderSeqEqual, num2str(uMaxVec(i)), num2str(threshVec(j)));
            fileSeqUnequal = sprintf("%sumax_%s_thresh_%s/", folderSeqUnequal, num2str(uMaxVec(i)), num2str(threshVec(j)));

            %population files
            populationFileJointEqual = sprintf("%sindividualGraphs/population_x2x3ratio0.5_concAconcBratio0.1_gammaEqual.mat", fileJointEqual);
            populationFileJointUnequal = sprintf("%sindividualGraphs/population_x2x3ratio0.5_concAconcBratio0.1_gammaUnequal.mat", fileJointUnequal);
            populationFileSeqEqual = sprintf("%sindividualGraphs/population_x2x3ratio0.5_gammaEqual.mat", fileSeqEqual);
            populationFileSeqUnequal = sprintf("%sindividualGraphs/population_x2x3ratio0.5_gammaUnequal.mat", fileSeqUnequal);
            
            %load in pop structs
            structJointEqual = load(populationFileJointEqual, "popStruct");
            structJointUnequal = load(populationFileJointUnequal, "popStruct");
            structSeqEqual = load(populationFileSeqEqual, "popStruct");
            structSeqUnequal = load(populationFileSeqUnequal, "popStruct");

            %establish best and worst for equal and unequal settings
            idealEqual = calculateIdealFromPopStruct(structJointEqual.popStruct);
            idealUnequal = calculateIdealFromPopStruct(structJointUnequal.popStruct);

            worstEqual = calculateWorstFromPopStruct(structJointEqual.popStruct);
            worstUnequal = calculateWorstFromPopStruct(structJointUnequal.popStruct);

            %heatmap file names
            heatmapFileJointEqual = sprintf("%ssimHeatmap.mat",fileJointEqual);
            heatmapFileJointUnequal = sprintf("%ssimHeatmap.mat",fileJointUnequal);
            heatmapFileSeqEqual = sprintf("%ssimHeatvector.mat",fileSeqEqual);
            heatmapFileSeqUnequal = sprintf("%ssimHeatvector.mat",fileSeqUnequal);
            
            %load in heatmaps
            heatmapJointEqual = load(heatmapFileJointEqual, "heatmapVals");
            heatmapJointUnequal = load(heatmapFileJointUnequal, "heatmapVals");
            heatmapSeqEqual = load(heatmapFileSeqEqual, "heatvectorVals");
            heatmapSeqUnequal = load(heatmapFileSeqUnequal, "heatvectorVals");

            jointVec = [max(heatmapJointEqual.heatmapVals, [], "all"), max(heatmapJointUnequal.heatmapVals, [], "all")];
            jointWorstVec = [min(nonzeros(heatmapJointEqual.heatmapVals), [], "all"), min(nonzeros(heatmapJointUnequal.heatmapVals), [], "all")];
            seqVec = [max(heatmapSeqEqual.heatvectorVals, [], "all"), max(heatmapSeqUnequal.heatvectorVals, [], "all")];
            seqWorstVec = [min(nonzeros(heatmapSeqEqual.heatvectorVals), [], "all"), min(nonzeros(heatmapSeqUnequal.heatvectorVals), [], "all")];

            figure((i-1)*length(threshVec)+j)
                hold on

                %best & worst lines
                line([0.95, 1.05], [idealEqual, idealEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([0.95, 1.05], [worstEqual, worstEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1, 1], [worstEqual, idealEqual], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1.95, 2.05], [idealUnequal, idealUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([1.95, 2.05], [worstUnequal, worstUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')
                line([2, 2], [worstUnequal, idealUnequal], 'LineWidth', 2.5, 'Color', '#bbb', 'HandleVisibility', 'off')

                %plot the various bits and bobs
                plot([0.975, 1.975], jointVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Best Joint')
                plot([0.9875, 1.9875], seqVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Best Sequential')
                plot([1.0125, 2.0125], jointWorstVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Worst Joint')
                plot([1.025, 2.025], seqWorstVec, 'LineStyle','none', 'Marker', '.', 'MarkerSize', 30, 'DisplayName', 'Worst Sequential')

                %axis settings
                maxIdeal = max(idealUnequal, idealEqual)+50;
                axis([0 3 0 maxIdeal])

                legend show

                ylabel("time (hours)", 'FontSize', 18)
                xticks([1 2])
                xticklabels({"Equal \gamma", "Unequal \gamma"});

                title(sprintf("uMax - %s Threshold - %s",num2str(uMaxVec(i)), num2str(threshVec(j))))

                set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
                set(gca, 'FontSize', 16)
                hold off

        end
    end
end

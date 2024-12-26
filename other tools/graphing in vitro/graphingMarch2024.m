normRed = 4827.996032;
normYellow = 38760.84921;
normGreen = 589752.078;
normBlue = 39207.43651;


worst_red = readtable("worstAnalysis.xlsx", sheet="red");
wRed_arr = table2array(worst_red);
worst_yellow = readtable("worstAnalysis.xlsx", sheet="yellow");
wYellow_arr = table2array(worst_yellow);
worst_green = readtable("worstAnalysis.xlsx", sheet="green");
wGreen_arr = table2array(worst_green);
worst_blue = readtable("worstAnalysis.xlsx", sheet="blue");
wBlue_arr = table2array(worst_blue);

time = wRed_arr(:, 1);
blueTime = wRed_arr(6:end, 1);
avg_wRed = wRed_arr(:, 3)/normRed;
avg_wYellow = wYellow_arr(:, 3)/normYellow;
avg_wGreen = wGreen_arr(:, 3)/normGreen;
avg_wBlue = wBlue_arr(:, 3)/normBlue;

figure(1) %worst figure
hold on;
plot(time, avg_wRed, 'LineWidth', 4, 'Color', '#A2142F');
plot(time, avg_wYellow, 'LineWidth', 4, 'Color', '#EDB120');
plot(time, avg_wGreen, 'LineWidth', 4, 'Color', '#77AC30');
plot(blueTime, avg_wBlue, 'LineWidth', 4, 'Color', '#0072BD');
xline([5, 10, 14, 18, 22], 'LineStyle','--', 'LineWidth', 2)
legend('No resistance', 'Afatinib resistance', 'Both resistance', 'Osimertinib resistance', 'Location','northwest')
set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 18);
title('Worst Regimen', 'FontSize', 24);
xlabel('time (days)', 'FontSize', 20);
ylabel('normalized fluorescent intensity', 'FontSize', 20);
hold off; 

afaStart_red = readtable("bestsAfaStartAnalysis.xlsx", sheet="red");
afaRed_arr = table2array(afaStart_red);
afaStart_yellow = readtable("bestsAfaStartAnalysis.xlsx", sheet="yellow");
afaYellow_arr = table2array(afaStart_yellow);
afaStart_green = readtable("bestsAfaStartAnalysis.xlsx", sheet="green");
afaGreen_arr = table2array(afaStart_green);
afaStart_blue = readtable("bestsAfaStartAnalysis.xlsx", sheet="blue");
afaBlue_arr = table2array(afaStart_blue);

avg_afaRed = afaRed_arr(:, 3)/normRed;
avg_afaYellow = afaYellow_arr(:, 3)/normYellow;
avg_afaGreen = afaGreen_arr(:, 3)/normGreen;
avg_afaBlue = afaBlue_arr(:, 3)/normBlue;

figure(2) %worst figure
hold on;
plot(time, avg_afaRed, 'LineWidth', 4, 'Color', '#A2142F');
plot(time, avg_afaYellow, 'LineWidth', 4, 'Color', '#EDB120');
plot(time, avg_afaGreen, 'LineWidth', 4, 'Color', '#77AC30');
plot(blueTime, avg_afaBlue, 'LineWidth', 4, 'Color', '#0072BD');
xline([5, 10, 14, 18, 22], 'LineStyle','--', 'LineWidth', 2)
legend('No resistance', 'Afatinib resistance', 'Both resistance', 'Osimertinib resistance', 'Location','northwest')
set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 18);
title('Afatinib Start Alt Regimen', 'FontSize', 24);
xlabel('time (days)', 'FontSize', 20);
ylabel('normalized fluorescent intensity', 'FontSize', 20);
hold off;

osiStart_red = readtable("bestsOsiStartAnalysis.xlsx", sheet="red");
osiRed_arr = table2array(osiStart_red);
osiStart_yellow = readtable("bestsOsiStartAnalysis.xlsx", sheet="yellow");
osiYellow_arr = table2array(osiStart_yellow);
osiStart_green = readtable("bestsOsiStartAnalysis.xlsx", sheet="green");
osiGreen_arr = table2array(osiStart_green);
osiStart_blue = readtable("bestsOsiStartAnalysis.xlsx", sheet="blue");
osiBlue_arr = table2array(osiStart_blue);

avg_osiRed = osiRed_arr(:, 3)/normRed;
avg_osiYellow = osiYellow_arr(:, 3)/normYellow;
avg_osiGreen = osiGreen_arr(:, 3)/normGreen;
avg_osiBlue = osiBlue_arr(:, 3)/normBlue;

figure(3) %worst figure
hold on;
plot(time, avg_osiRed, 'LineWidth', 4, 'Color', '#A2142F');
plot(time, avg_osiYellow, 'LineWidth', 4, 'Color', '#EDB120');
plot(time, avg_osiGreen, 'LineWidth', 4, 'Color', '#77AC30');
plot(blueTime, avg_osiBlue, 'LineWidth', 4, 'Color', '#0072BD');
xline([5, 10, 14, 18, 22], 'LineStyle','--', 'LineWidth', 2)
legend('No resistance', 'Afatinib resistance', 'Both resistance', 'Osimertinib resistance', 'Location','northwest')
set(gcf, 'units', 'points', 'position', [10, 10, 1080, 864]); %points = 1/72", ~8" here?
set(gca, 'FontSize', 18);
title('Osimertinib Start Alt Regimen', 'FontSize', 24);
xlabel('time (days)', 'FontSize', 20);
ylabel('normalized fluorescent intensity', 'FontSize', 20);
hold off;
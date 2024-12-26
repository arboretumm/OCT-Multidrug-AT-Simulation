function megarun_24_0611_uMaxkR()
   uMaxVec = [0.5, 0.75, 1, 10, 100];
   threshVec = [0.25, 0.5, 0.8, 0.99];

   for i=1:length(uMaxVec)
       for j=1:length(threshVec)
           disp("this is run")
           disp((i-1)*(length(uMaxVec))+j)
           multisimRun_24_0611_jointWithSwitchOldGamma_heatmapSave(uMaxVec(i), threshVec(j));
       end
   end
end
function ratio = compareResponseAtPositionBetweenTwoRastermaps(rmaps1,rmaps2,position)

if nargin<3
   position = 22; 
end

% Parameters
params.width = 3;

% Init
ratio = [];

% Create average tuning and divide into odd and even laps
avg_rmaps_on = sb.generate.averagedTuningFromRasterMaps(rmaps1);
avg_rmaps_off = sb.generate.averagedTuningFromRasterMaps(rmaps2);

% Compute ratio of the amplitude in ON laps vs OFF laps for putative cells
for c = 1:size(avg_rmaps_on,1)
    
   [~,peak_position] = max(avg_rmaps_on(c,position:position+20));
   peak_position = peak_position + position;
   
   amplitude_on = nanmean(avg_rmaps_on(c,peak_position-params.width:peak_position+params.width));
   amplitude_off = nanmean(avg_rmaps_off(c,peak_position-params.width:peak_position+params.width));

   ratio(c) = ((amplitude_off/amplitude_on)-1)*100;
   
   if isinf(ratio(c)) || ratio(c)>1000
       ratio(c)= NaN;
   end

end






end
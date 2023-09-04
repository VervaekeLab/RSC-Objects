%% Use the sessions where the video captured the interaction between landmarks and whisker to see whe
% Use the sessions where the video captured the interaction between 
% landmarks and whisker to see which bin would be the earliest where the
% mouse touches the landmark.

landmark_side = [1,0,1,1]; % 1 = left
figure(1);clf;

for s = 1:4
    
   % Get landmark signal
   if landmark_side(s)
       landmark = mData(s).sData.whiskers.left_motor_sig;
   else
       landmark = mData(s).sData.whiskers.right_motor_sig;
   end
      
   % Convert landmark signal
   landmark(landmark<70) = 1000;
   landmark(landmark<1000) = 0;
   
   landmark = smoothdata(landmark,'movmean',10);
   landmark(landmark>0) = 1;
   
   landmark = landmark(mData(s).sData.whiskers.led_onsets);
   
   % Get behavioral data and create rastermap
   pos = mData(s).sData.behavior.wheelPosDsBinned;
   laps = mData(s).sData.behavior.wheelLapDsBinned;
   
   if min(laps)<0
       laps = laps + abs(min(laps))+1;
   end
   
   rmap = sb.generate.rasterMaps(pos,laps,landmark);

   landmark_onsets = find(nanmean(rmap));
   
   subplot(2,2,s);
   imagesc(rmap);
 
   hold on;
   line([],[]),
   first_landmark = find(ismember(landmark_onsets,30:40));
   first_landmark = first_landmark(1);
   second_landmark = find(ismember(landmark_onsets,70:80));
   second_landmark = second_landmark(1);
   
   first_onsets(s) = first_landmark;
   second_onsets(s) = second_landmark;
   
end
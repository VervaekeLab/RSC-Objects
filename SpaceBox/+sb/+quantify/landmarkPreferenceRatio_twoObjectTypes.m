function [identityIdxA,identityIdxB,ratio] = landmarkPreferenceRatio_twoObjectTypes(rmaps,paramsA,paramsB)
% Compute the ratio of preference for each cell between the two landmark
% positions. If the first landmark is twice as big as the second, it will
% be given a score of -100%, indicating a 100% increase difference towards 
% where - means it is the first landmark. 


% Compute the ratio between the two, so basically dividing
% first and second landmark amplitude
ratio = [];


for c = 1:size(rmaps,3)
    
    % Create position tuning map for the cell and normalize it
    avg_tuning = normalize(nanmean(rmaps(:,:,c)),'range',[0,1]);

    
    % object A find both peak responses
    % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
    [~,first_peakA] = max(avg_tuning(paramsA.first_landmark_bins));
    [~,second_peakA] = max(avg_tuning(paramsA.second_landmark_bins));

    % Set correct position
    first_peakA = first_peakA + paramsA.first_landmark_bins(1);
    second_peakA = second_peakA + paramsA.second_landmark_bins(1);

    % Find the average response around the two peaks in a 10 cm window
    first_responseA = nanmean(avg_tuning(first_peakA-3:first_peakA+3));
    second_responseA = nanmean(avg_tuning(second_peakA-3:second_peakA+3));

    % Find biggest and smallest
%     [biggest_responseA,biggest_indA] = max([first_responseA,second_responseA]);
%     smallest_responseA = min([first_responseA,second_responseA]);
%     
%     
    
     % object B find both peak responses
    [~,first_peakB] = max(avg_tuning(paramsB.first_landmark_bins));
    [~,second_peakB] = max(avg_tuning(paramsB.second_landmark_bins));

    % Set correct position
    first_peakB = first_peakB + paramsB.first_landmark_bins(1);
    second_peakB = second_peakB + paramsB.second_landmark_bins(1);

    % Find the average response around the two peaks in a 10 cm window
    first_responseB = nanmean(avg_tuning(first_peakB-3:first_peakB+3));
    second_responseB = nanmean(avg_tuning(second_peakB-3:second_peakB+3));

    % Find biggest response of not prominent object and smallest of most prominent object 

    [~ ,biggest_ind] = max([first_responseA,second_responseA,first_responseB,second_responseB]);
    identityIdxA(c) = NaN;
    identityIdxB(c) = NaN;
    if biggest_ind < 3 % object A is dominant
        identityIdxA(c) = min([first_responseA,second_responseA])-max([first_responseB,second_responseB]);
        ratio(c) = (nanmean([first_responseB,second_responseB])/nanmean([first_responseA,second_responseA])-1)*100;
        
%         identityIdx(c) = - identityIdx(c);
        ratio(c) = -ratio(c);
    else    % object B is dominant
        identityIdxB(c) = min([first_responseB,second_responseB])-max([first_responseA,second_responseA]);
        ratio(c) = (nanmean([first_responseB,second_responseB])/nanmean([first_responseA,second_responseA])-1)*100;
    end
    
    
    
    
    
    
    
    
    
    
% 
%     ratio(c) = ((biggest_response/smallest_response)-1) * 100;
%     
%     % If the first peak is the biggest, then set ratio to negative.
%     if biggest_ind == 1
%         ratio(c) = - ratio(c);
%     end
    
       
end

end
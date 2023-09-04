function [identityIdxA,identityIdxB,ratio] = landmarkPreferenceRatio_twoObjectTypesShuffle(rmaps,paramsA,paramsB)
% Compute the ratio of preference for each cell between the two landmark
% positions. If the first landmark is twice as big as the second, it will
% be given a score of -100%, indicating a 100% increase difference towards 
% where - means it is the first landmark. 


% Compute the ratio between the two, so basically dividing
% first and second landmark amplitude
ratio = [];


for c = 1:size(rmaps,3)
    
    rmap = rmaps(:,:,c);
   
    shuffled_avg = zeros(1000,size(rmap,2));
    for i = 1:1000

        % Init zeroed matrix
        shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

        % Shuffle for each lap
        for l = 1:size(rmap,1)
           shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(rmap,2))); 
        end

         % Create position tuning map for the cell and normalize it

        shuffled_avg(i,:) = normalize(nanmean(shuffled_rmap),'range',[0,1]);         

        
    
    % object A find both peak responses
    % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
    [~,first_peakA] = max(shuffled_avg(i,paramsA.first_landmark_bins));
    [~,second_peakA] = max(shuffled_avg(i,paramsA.second_landmark_bins));

    % Set correct position
    first_peakA = first_peakA + paramsA.first_landmark_bins(1);
    second_peakA = second_peakA + paramsA.second_landmark_bins(1);

    % Find the average response around the two peaks in a 10 cm window
    first_responseA = nanmean(shuffled_avg(i,first_peakA-3:first_peakA+3));
    second_responseA = nanmean(shuffled_avg(i,second_peakA-3:second_peakA+3));

    
     % object B find both peak responses
    [~,first_peakB] = max(shuffled_avg(i,paramsB.first_landmark_bins));
    [~,second_peakB] = max(shuffled_avg(i,paramsB.second_landmark_bins));

    % Set correct position
    first_peakB = first_peakB + paramsB.first_landmark_bins(1);
    second_peakB = second_peakB + paramsB.second_landmark_bins(1);

    % Find the average response around the two peaks in a 10 cm window
    first_responseB = nanmean(shuffled_avg(i,first_peakB-3:first_peakB+3));
    second_responseB = nanmean(shuffled_avg(i,second_peakB-3:second_peakB+3));

    % Find biggest response of not prominent object and smallest of most prominent object 

    [~ ,biggest_ind] = max([first_responseA,second_responseA,first_responseB,second_responseB]);
    identityIdxA(i,c) = NaN;
    identityIdxB(i,c) = NaN;
    if biggest_ind < 3 % object A is dominant
        identityIdx(i,c) = min([first_responseA,second_responseA])-max([first_responseB,second_responseB]);
        ratio(i,c) = (nanmean([first_responseB,second_responseB])/nanmean([first_responseA,second_responseA])-1)*100;
        
        identityIdxA(i,c) = identityIdx(i,c);
        ratio(i,c) = -ratio(i,c);
    else    % object B is dominant
        identityIdxB(i,c) = min([first_responseB,second_responseB])-max([first_responseA,second_responseA]);
        ratio(i,c) = (nanmean([first_responseB,second_responseB])/nanmean([first_responseA,second_responseA])-1)*100;
    end
    
    end
    
    
    
       
end

end
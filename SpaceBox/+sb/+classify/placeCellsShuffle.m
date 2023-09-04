function [pcs] = placeCellsShuffle(rmaps,varargin)

% set seed to make reproducible
rng(5)

pcs = [];
for c = 1:size(rmaps,3)
     rmap = normalize(rmaps(:,:,c),'range',[0 1]);
    response = nanmean(rmap);

    shuffled_avg = zeros(1000,size(response,2));
    for i = 1:1000

        % Init zeroed matrix
        shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

        % Shuffle for each lap
        for l = 1:size(rmap,1)
            shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(response,2)));
        end

        shuffled_avg(i,:) = nanmean(shuffled_rmap);

    end


    baseline_threshold = prctile(shuffled_avg,99);
    
    significant_pos = find(response>baseline_threshold);
    switch_peak  = [find(diff(significant_pos)>1),length(significant_pos)];%[find(diff(significant_pos)> 4),length(significant_pos)];
    length_peaks = [switch_peak(1),diff(switch_peak)];
    
    % correct for circularity of treadmill
    if ~isempty(significant_pos)
        if significant_pos(1) == 1 && significant_pos(end) == 105
           length_peaks(1) = length_peaks(1)+length_peaks(end);
           length_peaks(end)= [];
        end
    end
    
    n_peaks(c)   = length(find(length_peaks>3));
    
   
    if ~isempty(find(n_peaks(c) > 0))
        pcs = [pcs,c];
    end

    
    
end
    
pcs = unique(pcs);


end



function [lcs] = landmarkCellsBasedOnTemplateMatching(rmaps,varargin)

%% Parameters
params.visualize = false;

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Find landmark cells based on rmaps
% Init
lcs = [];

% Create a template
template_cell = zeros(1,105);
template_cell([32,72]) = 1;
template_cell = normalize(smoothdata(template_cell,'gaussian',10),'range',[0,1]);

for r = 1:size(rmaps,3)
    
    % Create lap-averaged map of roi
    rmap = rmaps(:,:,r);
    avg_map = nanmean(rmap);
    
    % Smooth the averaged map
    avg_map_smoothed = normalize(smoothdata(avg_map,'movmean',20),'range',[0,1]);
    
    correlations = [];
    
    for shift = 1:105
       
        template_shifted = circshift(template_cell,shift);
        corr = corrcoef(avg_map_smoothed,template_shifted);
        correlations(shift) = corr(1,2);
        
    end
    
    % Find the optimal shift
    [optimal_correlation,optimal_shift] = (max(correlations));
    template_shifted = circshift(template_cell,optimal_shift);
    
    % Shuffling
    response = avg_map;

    % Shuffle 1000 times
    shuffled_avg = zeros(1000,1);
    for i = 1:1000

        % Init zeroed matrix
        shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

        % Shuffle for each lap
        for l = 1:size(rmap,1)
           shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(response,2))); 
        end
        
        corr = corrcoef(nanmean(shuffled_rmap),template_shifted);

        shuffled_avg(i,1) = corr(1,2);            

    end

    baseline_correlation = prctile(shuffled_avg,99);
    
    if optimal_correlation > baseline_correlation
        lcs = [lcs,r];
    end
    
    
end

end
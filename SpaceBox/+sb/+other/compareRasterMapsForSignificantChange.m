function significant_response = compareRasterMapsForSignificantChange(rmap1,rmap2, varargin)
%% compareRasterMapsForSignificantChange
% Compare rmap1 and rmap2 to find out if there is a significant increase in
% activity in rmap2 (response) compared to rmap1 (baseline)
% INPUT
%   -- Required
%   rmap1:
%   rmap2:
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   significant_response: Boolean


%% Default parameters
params = struct();
params.method = 'shuffle';

% Print all variables if "help" is only input 
try
    if inp == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

% Find the mean and std in baseline
baseline_mean = nanmean(rmap1(:));
baseline_std = nanstd(rmap1(:));
threshold = 3; % number of standard deviations above mean for a significant response to be detected

significant_response = zeros(1,size(rmap2,1));

switch params.method
    case 'mean'
        
        baseline = zeros(1,size(rmap1,1));
        for l = 1:size(rmap1,1)
            baseline(l) = nanmean(rmap1(l,:));
        end
        
        threshold = nanmean(baseline) + (2*nanstd(baseline));
        
        response = nanmean(nanmean(rmap2(:,20:51)));
        
        if response > threshold
           significant_response = 1;
        else
            significant_response = 0;
        end
       
        
    case 'max'
        for l = 1:size(rmap2,1)
            if max(smoothdata(rmap2(l,:),'gaussian',3)) > (baseline_mean + baseline_std*threshold)
                significant_response(l) = 1;
            end
        end
        
    case 'signrank'
        
        
        response = nanmean(rmap2(:,20:40),2);
        baseline = nanmean(rmap1(:,end-62:end),2);
        
        result = signrank(baseline,response,'tail','left');
        if result < 0.01
            significant_response = 1;
        else
            significant_response = 0;
        end
        
        
    case 'shuffle'
        
        rmap = [rmap1,rmap2];
        response = nanmean(rmap);
        
        % Shuffle 1000 times
        shuffled_avg = zeros(1000,155);
        for i = 1:1000
            
            % Init zeroed matrix
            shuffled_rmap = zeros(size(rmap,1),size(rmap,2));
            
            % Shuffle for each lap
            for l = 1:size(rmap,1)
               shuffled_rmap(l,:) = circshift(rmap(l,:),randi(155)); 
            end
            
            shuffled_avg(i,:) = nanmean(shuffled_rmap);            
            
        end
        
        baseline_threshold = prctile(shuffled_avg,95);
       
        if nanmean(response(110:130)) > nanmean(baseline_threshold(110:130))
            significant_response = 1;
        else
            significant_response = 0;
        end
        
end

% switch params.method
%     case 'mean'
%         
%         % Reliability
%         reliability = 0.25;
% 
%         if sum(significant_response)/length(significant_response) > reliability
%             significant_response = 1;
%         else
%             significant_response = 0;
%         end
% 
% end

end









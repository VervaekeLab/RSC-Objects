function reliability = touchReliability(rmaps,varargin)
%% touchReliability
% Quantify the reliability at landmark positions from a rastermap.
%
% INPUT
%   -- Required
%   rmaps: NxMxK matrix with K rastermaps of cells, N is lap, M is position
%   bin.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   reliability: The reliability score for each rastermap.
%

%% Default parameters
params = struct();
params.reliability_bins = {[23:50],[63:90]}; % What bins to be looking for activity
params.reliability_method = "zscoring";
params.subtract_baseline = false;

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


% For each rastermap compute reliabilty
reliability = zeros(size(rmaps,3),length(params.reliability_bins)+1);

for c = 1:size(rmaps,3)
    
    
    switch params.reliability_method
        
        case "zscoring"
           
            % Compute zscore normalized map
            z_value = 2; % Minimum z score value to be considered significant
            normalized_rmap = normalize(rmaps(:,:,c),2,'zscore');
            normalized_rmap(normalized_rmap<z_value) = 0;
            
%             % Plot cell
%             figure(1);
%             clf;
%             subplot(1,3,1);
%             imagesc(rmaps(:,:,c));
%             subplot(1,3,2);
%             imagesc(normalize(rmaps(:,:,c),2,'zscore'),[2,5]);
%             subplot(1,3,3);
%             imagesc(normalized_rmap,[0,z_value]);
            
            % Compute reliability for each landmark
            
            for l = 1:length(params.reliability_bins)
               
                bins = params.reliability_bins{l};
                
                % Extract only the  part of the map where the landmark is
                current_map = normalized_rmap(:,bins);
                
                % Count each laps with significant response at landmark
                significant_laps = 0;
                for lap = 1:size(current_map,1)
                    if nanmean(current_map(lap,:)) > 0
                       significant_laps = significant_laps + 1; 
                    end
                end
                
                reliability(c,l) = significant_laps / size(current_map,1);
                
            end
            
            % Compute a baseline reliability
            current_map = normalized_rmap(:,5:23);
            significant_laps = 0;
            
            for lap = 1:size(current_map,1)
                if nanmean(current_map(lap,:)) > 0
                    significant_laps = significant_laps + 1;
                end
            end
            
            reliability(c,l+1) = significant_laps / size(current_map,1);
            %reliability(c,:)
            %fprintf('Corrected: %.2f and %.2f\n',reliability(c,1)-reliability(c,3),reliability(c,2)-reliability(c,3));
    end
    
    
    
end

% Subtract baseline?
if params.subtract_baseline
   reliability(:,1) = reliability(:,1) - reliability(:,end);
   reliability(:,2) = reliability(:,2) - reliability(:,end);
else
    reliability(:,end) = [];
end


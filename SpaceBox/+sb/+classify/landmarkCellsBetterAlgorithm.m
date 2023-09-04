function [lcs, not_lcs] = landmarkCellsBetterAlgorithm(rmaps, varargin)
%% landmarkCellsBetterAlgorithm
% Improved method for classifying cells as landmark cells
% INPUT
%   -- Required
%
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%
%
% Written May 2021


%% Default parameters
params = struct();
params.method = 'amplitude_deconv';

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


%% Detection

lcs = [];

switch params.method
    
    % Pick cells that have an amplitude larger than some baseline at both
    % landmark positions
    
    case 'amplitude_deconv'
        
        first_landmark = 23:43;
        second_landmark = 63:83;
        
        % For each cell
        for c = 1:size(rmaps,3)
            
            % Get the cells rastermap
            rmap = rmaps(:,:,c);
            average_tuning = nanmean(rmap,1);
            
            [~,sort_ind] = sort(average_tuning);
            
            baseline_rmap = rmap(:,sort_ind(1:10));
            baseline_mean = nanmean(baseline_rmap(:));
            baseline_std = nanstd(baseline_rmap(:));
            
            threshold = baseline_mean + (baseline_std*4);
            
            % For each lap compute reliability
            reliability_first = [];
            reliability_second = [];
            
            for l = 1:size(rmap,1)
                
                % First landmark
                if max(rmap(l,first_landmark)) > threshold
                    reliability_first(l) = 1;
                else
                    reliability_first(l) = 0;
                end
                
                % Second landmark
                if max(rmap(l,second_landmark)) > threshold
                    reliability_second(l) = 1;
                else
                    reliability_second(l) = 0;
                end
                
                
            end
            
            % Compute reliability score for both landmarks for the cell
            reliability_first = sum(reliability_first) / length(reliability_first);
            reliability_second = sum(reliability_second) / length(reliability_second);
            
            
            % If the cell has a reliability score above 20% for both
            % landmarks, it is considered a landmark cell.

            if reliability_first > 0.33 && reliability_second > 0.33
               lcs(c) = 1;
            else
                lcs(c) = 0;
            end
    
        end
        
    case 'amplitude_dff'
        
        first_landmark = 20:60;
        second_landmark = 60:90;
        
        % For each cell
        for c = 1:size(rmaps,3)
            
            % Get the cells rastermap
            rmap = rmaps(:,:,c);
            average_tuning = nanmean(rmap,1);
            
            [~,sort_ind] = sort(average_tuning);
            
            baseline_rmap = rmap(:,sort_ind(1:10));
            baseline_mean = nanmean(baseline_rmap(:));
            baseline_std = nanstd(baseline_rmap(:));
            
            threshold = baseline_mean + (baseline_std*3);

            % For each lap compute reliability
            reliability_first = [];
            reliability_second = [];
            
            for l = 1:size(rmap,1)
                
                % First landmark
                if max(rmap(l,first_landmark)) > threshold
                    reliability_first(l) = 1;
                else
                    reliability_first(l) = 0;
                end
                
                % Second landmark
                if max(rmap(l,second_landmark)) > threshold
                    reliability_second(l) = 1;
                else
                    reliability_second(l) = 0;
                end
                
                
            end
            
            % Compute reliability score for both landmarks for the cell
            reliability_first = sum(reliability_first) / length(reliability_first);
            reliability_second = sum(reliability_second) / length(reliability_second);
            
            
            % If the cell has a reliability score above 20% for both
            % landmarks, it is considered a landmark cell.

            if reliability_first > 0.33 && reliability_second > 0.33
               lcs(c) = 1;
            else
                lcs(c) = 0;
            end
            
            
        end
        
        
        % Look at the correlation between a cells average tuning and the
        % optimal target tuning
    case 'correlation'
        
        % Target cell response
        target = zeros(1,105);
        target(36) = 1;
        target(76) = 1;
        
        target = smoothdata(target,'gaussian',20);
        target = normalize(target,'range',[0,1]);
        
        % For each cells
        corr = [];
        for c = 1:size(rmaps,3)
            
            % Get average response for the cell
            rmap = rmaps(:,:,c);
            average_tuning = nanmean(rmap,1);
            
            % Smooth average
            average_tuning = smoothdata(average_tuning,'gaussian',10);
            response = normalize(average_tuning,'range',[0,1]);
            
            % Compute correlation
            correlation = corrcoef(target,response);
            corr(c) = correlation(2,1);
        end
        
        case 'amplitude_dff_new'
        
        first_landmark = 20:60;
        second_landmark = 60:90;
        
        % For each cell
        for c = 1:size(rmaps,3)
            
            % Get the cells rastermap
            rmap = rmaps(:,:,c);
            average_tuning = nanmean(rmap,1);
            
            [~,sort_ind] = sort(average_tuning);
            
            baseline_rmap = rmap(:,sort_ind(1:6));
            baseline_mean = nanmean(baseline_rmap(:));
            baseline_std = nanstd(baseline_rmap(:));
            
            threshold = baseline_mean + (baseline_std*3);
            
            average_tuning_mod = average_tuning;
            average_tuning_mod(average_tuning_mod<threshold) = 0;
            [~,peak_position] = findpeaks(average_tuning_mod,'NPeaks',2,'MinPeakDist',20);
            
            lcs(c) = 0;
            if length(peak_position) == 2
                if sum(ismember(first_landmark,peak_position))
                    if sum(ismember(second_landmark,peak_position))
                        
                        lcs(c) = 1;
                        
                    end
                end
            end
            
%             % For each lap compute reliability
%             reliability_first = [];
%             reliability_second = [];
%             
%             for l = 1:size(rmap,1)
%                 
%                 % First landmark
%                 if max(rmap(l,first_landmark)) > threshold
%                     reliability_first(l) = 1;
%                 else
%                     reliability_first(l) = 0;
%                 end
%                 
%                 % Second landmark
%                 if max(rmap(l,second_landmark)) > threshold
%                     reliability_second(l) = 1;
%                 else
%                     reliability_second(l) = 0;
%                 end
%                 
%                 
%             end
%             
%             % Compute reliability score for both landmarks for the cell
%             reliability_first = sum(reliability_first) / length(reliability_first);
%             reliability_second = sum(reliability_second) / length(reliability_second);
%             
%             
%             % If the cell has a reliability score above 20% for both
%             % landmarks, it is considered a landmark cell.
% 
%             if reliability_first > 0.33 && reliability_second > 0.33
%                lcs(c) = 1;
%             else
%                 lcs(c) = 0;
%             end
            
            
        end
             
end


% Get correct index number
not_lcs = find([lcs==0]);
lcs = find(lcs);


end

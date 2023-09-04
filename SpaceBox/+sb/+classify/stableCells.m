function rmaps = stableCells(rmaps,varargin)


% NOT FINISHED


%% stableCells
% Exclude non stable cells based on the rastermap given as input and
% parameters set in the function. So, fex, you can find landmark cells that
% have a high stability score to do further analysis with.
% INPUT
%   -- Required
%   rmaps: rastermap for all cells, with x,y,z being lap, pos, cell.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   rmaps: rastermap for all cells that has passed the criteria.
%



%% Default parameters
params = struct();
params.criteria = "standard"; % The stability criteria that sets the method to find cells.
params.peaks = [36,76]; % Positions where peaks must be present in rastermap
params.peak_width = 5; % The acceptable +/- difference from the peaks to include 

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


cors = [];
for c = 1:size(rmaps,3)
    
    current_map = rmaps(:,:,c);
    %sb.plot.rasterMaps(current_map);
    
    
    % Autocorrelation between gaussian and signal?
    average = nanmean(current_map);
    
    signal = zeros(1,105);
    signal(36) = 1;
    %signal(76) = 1;
    signal = signal(25:95);
    signal = smoothdata(signal,'gaussian',20);
    average = average(25:95);
    
    cors(c) = max(xcorr(average,signal));
    
    
%     figure(29);
%     clf;
%     plot(normalize(average,'range',[0,1]));
%     hold on
%     plot(normalize(signal,'range',[0,1]))
%     plot(xcorr(average,signal))
%     
%     % Look at each peak
%     for p = params.peaks
%         poi = params.peaks(p)-params.peak_widthparams.peaks(p)+params.peak_width;
%         
%         peak_map = current_map(:,poi);
%         
%     end
    
    
   
    
end

[~,sort_ind] = sort(cors);
sb.plot.rasterMaps(rmaps(:,:,sort_ind),'rmap_title',num2cell(cors(sort_ind)));


end



function lcs = landmarkCellsStrict(rmaps,varargin)

% NOT FINISHED


%% landmarkCellsStrict
% Conservative method for classifying rastermaps as landmark cells.
% INPUT
%   -- Required
%   rmaps: rastermap of responses.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   lcs: list of inds of the rmaps which are deemed landmark cells.
%


%% Default parameters
params = struct();

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

lcs = [];
for c = 1:size(rmaps,3)
    roi = rmaps(:,:,c);
    roi_baseline = roi(:,5:30);
    baseline_mean = nanmean(roi_baseline(:));
    baseline_std = nanstd(roi_baseline(:));

    roi_response = roi(:,31:40);

    response_mean = nanmean(roi_response(:));
    response_std = nanstd(roi_response(:));


    if response_mean > baseline_mean + (3*baseline_std)
        %sb.plot.rasterMaps(roi);
        lcs = [lcs,c];
    end
end



% Try to find the average rastermaps that is the most correlated with a
% double peaked gaussian
rmaps_smoothed = smoothdata(rmaps,'gaussian',3);
target_signal = zeros(1,105);
target_signal(36) = 1;
target_signal(76) = 1;
target_signal = normalize(smoothdata(target_signal,'gaussian',30),'range',[0,1]);

average_map = normalize(squeeze(nanmean(rmaps_smoothed,1)),'range',[0,1]);
correlations = [];

for c = 1:size(rmaps,3)
    corrs = corrcoef(average_map(:,c),target_signal);
    
    correlations = [correlations, corrs(2,1)];
end

% 
% figure(22);
% clf;
% subplot(2,1,1);
% average_map = normalize(squeeze(nanmean(rmaps_smoothed,1)),'range',[0,1]);
% c = 37;
% plot(average_map(:,c));
% hold on
% plot(target_signal)
% subplot(2,1,2);
% plot(xcorr(average_map(:,c),target_signal))
% 
% 
% sb.plot.rasterMaps(rmaps(:,:,find([correlations>0.5])))



% Compare correlation of odd vs even laps
average_map_odd = normalize(squeeze(nanmean(rmaps(1:2:end,:,:),1)),'range',[0,1]);
average_map_even = normalize(squeeze(nanmean(rmaps(2:2:end,:,:),1)),'range',[0,1]);
correlations = [];

for c = 1:size(rmaps,3)
    corrs = corrcoef(average_map_odd(:,c),average_map_even(:,c));
    correlations = [correlations, corrs(2,1)];
end

figure(22);clf;plot(average_map_odd(:,1)); hold on; plot(average_map_even(:,1))

figure(23);
clf;
histogram(correlations)




end


function [mData,data] = rsc_jitter()
  %% analysis of tactile object jittering
%
% INPUT
%   -- Required: dirct
%
% OUTPUT
%      mData containing rmaps
%      data with classification
%

% please add loading direction to sData structs here!
dirct = '';
% Recordings to include
data.sessionIDs = { 'm1551-20230308-03', ... # Landmarks removed after 41 laps. Ran 96 laps total
                    'm1551-20230309-01', ... # Landmarks removed after 41 laps. Ran 136 laps total
                    'm1553-20230308-01', ... # Landmarks removed after 30 laps. Ran 44 laps total
                    'm1553-20230309-01', ...
                    'm1554-20230308-01',...
                    'm1554-20230309-01'
                    };
                
                
% % Load all sessions

for f = 1:length(data.sessionIDs)
    
    load(fullfile(dirct,data.sessionIDs{f},'sData.mat'));
    mData(f).sData = sData;
end



% Initiate some variables
% data.n_cells_pos1 = [];
data.n_lcs_pos1 = 0;
data.n_lcs_pos2 = 0;
data.n_lcs_pos3 = 0;
data.n_lcs_pos4 = 0;

max_z_value = 3;

% Create rastermaps for each session
for s = 1:length(mData)
    
    mData(s).sData.imdata.roiSignals(2).('dffSubtractedNpil') = mData(s).sData.imdata.roiSignals(2).dff;
    
    mData(s).sData.rmaps = sb.classify.split_rmaps_jitter(mData(s).sData);

    % identify landmark and place cells in on trials 
 
    params.first_landmark_bins = 15:29; params.second_landmark_bins = 48:62;
    mData(s).rmaps_first = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_firstpair,params);  
    params.first_landmark_bins = 25:39; params.second_landmark_bins = 58:72;    
    mData(s).rmaps_second = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_secondpair,params); 
    params.first_landmark_bins = 35:48; params.second_landmark_bins = 68:82;   
    mData(s).rmaps_third = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_thirdpair,params); 
    params.first_landmark_bins = 45:59; params.second_landmark_bins = 78:92;  
    mData(s).rmaps_fourth = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_fourthpair,params);  

   data.n_lcs_pos1 = data.n_lcs_pos1 +length(mData(s).rmaps_first);
   data.n_lcs_pos2 = data.n_lcs_pos2 + length(mData(s).rmaps_second);
   data.n_lcs_pos3 = data.n_lcs_pos3 + length(mData(s).rmaps_third);
   data.n_lcs_pos4 = data.n_lcs_pos4 + length(mData(s).rmaps_fourth);
    
   lcs_in_12 = intersect(mData(s).rmaps_first,mData(s).rmaps_second);
   lcs_in_123 =  intersect(lcs_in_12,mData(s).rmaps_third);
   data.lcs_in_all{s} =  intersect(lcs_in_123,mData(s).rmaps_fourth);
   data.n_lcs_in_all(s) = length(data.lcs_in_all{s});
   data.prct_lcs_in_all(s) = length(data.lcs_in_all{s})/size(mData(s).sData.rmaps.deconv_firstpair,3);
   

end
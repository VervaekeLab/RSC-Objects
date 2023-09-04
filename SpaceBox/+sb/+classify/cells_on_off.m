function rmaps = cells_on_off(sData,varargin)
%% Cells
% Find all spatially tuned cells in a sData struct.
% INPUT
%   -- Required
%   sData: Struct of session data.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   rmaps: Struct containing rastermaps for all cell types, and a subfield
%       called "ind" which contains the index of each cell class.

%% Default parameters
params = struct();
params.manipulated_motor = 1;
params.dff_signal = "dffSubtractedNpil"; % The name of the roiSignal field that contains the dff signal to be used. Either "dff" or "dffSubtractedNpil".

% Turn off displaying output
global verbose; verbose = false;
warning('off','all')

% Print all variables if "help" is only input 
try
    if sData == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Create rastermaps
rmaps = sb.generate.rasterMapsForSession(sData,"dff_signal", params.dff_signal,'manipulated_motor',params.manipulated_motor);

% Number of rois
n_rois = size(rmaps.dff,3);

%% Classify cell types
% Find active cells
acs = sb.classify.activeRoisAPT(sData.imdata.roiSignals(2).(params.dff_signal));

% % Find place cells

pcs = sb.classify.placeCellsShuffle(rmaps.deconv_motor_on);
% % % Find landmark cells and trace cells
lcs = sb.classify.landmarkCellsByShuffling(rmaps.deconv_motor_on);
tcs = sb.classify.landmarkCellsByShuffling(rmaps.deconv_motor_off);

pcs_off = sb.classify.placeCellsShuffle(rmaps.deconv_motor_off);


lcs = setdiff(lcs,tcs);


% Make place cells mutually exclusive from lcs and tcs
pcs = pcs(~ismember(pcs,lcs));
pcs = pcs(~ismember(pcs,tcs));
pcs_off = pcs_off(~ismember(pcs_off,tcs));

% Include only cells that show dFF activity
lcs = lcs(ismember(lcs,acs));
tcs = tcs(ismember(tcs,acs));
pcs = pcs(ismember(pcs,acs));
pcs_off=pcs_off(ismember(pcs_off,acs));
% Find all other cells that are not tuned
ocs = ones(1,n_rois);
ocs(unique([lcs,pcs,tcs])) = 0;
ocs = find(ocs);

% Landmark cells
rmaps.lcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,lcs);
rmaps.lcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,lcs);
rmaps.lcs.dff_motor_on = rmaps.dff_motor_on(:,:,lcs);
rmaps.lcs.dff_motor_off = rmaps.dff_motor_off(:,:,lcs);
rmaps.lcs.ind = lcs;

% Trace cells
rmaps.tcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,tcs);
rmaps.tcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,tcs);
rmaps.tcs.dff_motor_on = rmaps.dff_motor_on(:,:,tcs);
rmaps.tcs.dff_motor_off = rmaps.dff_motor_off(:,:,tcs);
rmaps.tcs.ind = tcs;

% Place cells
rmaps.pcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,pcs);
rmaps.pcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,pcs);
rmaps.pcs.dff_motor_on = rmaps.dff_motor_on(:,:,pcs);
rmaps.pcs.dff_motor_off = rmaps.dff_motor_off(:,:,pcs);
rmaps.pcs.ind = pcs;
rmaps.pcs_off.ind = pcs_off;
% Other cells
rmaps.ocs.deconv_motor_on = rmaps.deconv_motor_on(:,:,ocs);
rmaps.ocs.deconv_motor_off = rmaps.deconv_motor_off(:,:,ocs);
rmaps.ocs.dff_motor_on = rmaps.dff_motor_on(:,:,ocs);
rmaps.ocs.dff_motor_off = rmaps.dff_motor_off(:,:,ocs);
rmaps.ocs.ind = ocs;

% All rastermaps combined
rmaps.all.dff_motor_on = cat(3,rmaps.lcs.dff_motor_on,rmaps.tcs.dff_motor_on,rmaps.pcs.dff_motor_on);
rmaps.all.dff_motor_off = cat(3,rmaps.lcs.dff_motor_off,rmaps.tcs.dff_motor_off,rmaps.pcs.dff_motor_off);
rmaps.all.deconv_motor_on = cat(3,rmaps.lcs.deconv_motor_on,rmaps.tcs.deconv_motor_on,rmaps.pcs.deconv_motor_on);
rmaps.all.deconv_motor_off = cat(3,rmaps.lcs.deconv_motor_off,rmaps.tcs.deconv_motor_off,rmaps.pcs.deconv_motor_off);

% Add some metadata
rmaps.meta.n_cells = n_rois;
rmaps.meta.n_acs = length(acs);
rmaps.meta.n_pcs = length(pcs);
rmaps.meta.n_lcs = length(lcs);
rmaps.meta.n_tcs = length(tcs);
rmaps.meta.n_ocs = length(ocs);
rmaps.meta.n_pcs_off = length(pcs_off);

%fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nPure landmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\nAll landmark cells : %i (%.1f%%)\n',rmaps.meta.n_cells,rmaps.meta.n_acs,sum([rmaps.meta.n_pcs,rmaps.meta.n_lcs,rmaps.meta.n_tcs]),(sum([rmaps.meta.n_pcs,rmaps.meta.n_lcs,rmaps.meta.n_tcs])/rmaps.meta.n_acs)*100,rmaps.meta.n_pcs,(rmaps.meta.n_pcs/rmaps.meta.n_acs)*100,rmaps.meta.n_lcs,(rmaps.meta.n_lcs/rmaps.meta.n_acs)*100,rmaps.meta.n_tcs,(rmaps.meta.n_tcs/rmaps.meta.n_acs)*100,sum([rmaps.meta.n_lcs,rmaps.meta.n_tcs]),(sum([rmaps.meta.n_lcs,rmaps.meta.n_tcs])/rmaps.meta.n_acs)*100);
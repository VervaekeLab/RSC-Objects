function rmaps_ = classifyCells(sData,rmap_deconv)
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

% Number of rois
rng(2)
n_rois = size(rmap_deconv,3);

%% Classify cell types
% Find active cells
acs = sb.classify.activeRoisAPT(sData.imdata.roiSignals(2).('dffSubtractedNpil'));

% Find place cells
pcs = sb.classify.placeCellsShuffle(rmap_deconv);
lcs = sb.classify.landmarkCellsByShuffling(rmap_deconv(:,:,pcs));
lcs = pcs(lcs);
pcs = setdiff(pcs,lcs);

% Include only cells that show dFF activity
% lcs = lcs(ismember(lcs,acs));
% pcs = pcs(ismember(pcs,acs));

% Find all other cells that are not tuned
ocs = ones(1,n_rois);
ocs(unique([lcs,pcs])) = 0;
ocs = find(ocs);

% Landmark cells
rmaps_.lcs.deconv = rmap_deconv(:,:,lcs);
rmaps_.lcs.ind    = lcs;

% Place cells
rmaps_.pcs.deconv = rmap_deconv(:,:,pcs);
rmaps_.pcs.ind    = pcs;

% Other cells
rmaps_.ocs.deconv = rmap_deconv(:,:,ocs);
rmaps_.ocs.ind    = ocs;

% Add some metadata
rmaps_.meta.n_cells = n_rois;
rmaps_.meta.n_acs = length(acs);
rmaps_.meta.n_pcs = length(pcs);
rmaps_.meta.n_lcs = length(lcs);
rmaps_.meta.n_ocs = length(ocs);

function nmap = rasterMapByLandmarkOnset(rmap,onsets)
%% rasterMapByLandmarkOnset
%
% INPUT
%   -- Required
%   rmap: Rastermap matrix.
%   onsets: Array of [0, 1] of same length as dim 1 on rmap.
%
% OUTPUT
%   rmap: Same as input but now sorted
%
% Written by Andreas S Lande 2019


% Find out if rmap has 3 dims
n_maps = 1;

onsets = onsets(1:size(rmap,1)); % cut out exess data

% on and off laps
ons = find(onsets == 1);
offs = find(onsets == 0);

for m = 1:size(rmap,3)
    nmap(:,:,m) = [rmap(ons,:,m); rmap(offs,:,m)];
end

    
end
    
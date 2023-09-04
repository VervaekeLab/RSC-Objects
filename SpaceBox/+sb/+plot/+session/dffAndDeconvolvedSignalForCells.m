function dffAndDeconvolvedSignalForCells(sData,inds,varargin)
%% dffAndDeconvolvedSignalForCells
% Plot the deconvolved and dff signal for the list of cells given as input
% INPUT
%   -- Required
%       sData
%       inds: list of indexed number of all cells to be plotted
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%
%
% Written by Andreas S Lande 2020


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


%% Find the cells
signals_dff = [];
signals_deconv = [];

for c = 1:length(inds)
    signals_dff = [signals_dff; sData.imdata.roiSignals(2).dffSubtractedNpil(inds(c),:)];
    signals_deconv = [signals_dff; sData.imdata.roiSignals(2).deconv(inds(c),:)];
end


figure(2121);
clf;

plot(signals_dff,'LineWidth',2,'Color','k');







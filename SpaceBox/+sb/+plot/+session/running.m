function rmap = running(sData,varargin)
%% Running
% Plot the running profile of sData session.
% INPUT
%   -- Required
%       sData. 
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
%
% Written by Andreas S Lande 2020


%% Default parameters
params = struct();
params.createFigure = true;

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

%% Generate data
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
speed = sData.behavior.runSpeedDs;

% Plot
rmap = sb.generate.rasterMaps(pos,lap,speed,'smoothing',3);
if params.createFigure
    tit = sprintf('Running speed - %s',strrep(sData.sessionInfo.sessionID(1:17),'_','-'));
    sb.plot.rasterMaps(rmap, 'cmap','gray','rmap_title',{tit},'font_size',16, 'xlabel','Position (cm)','ylabel','Lap #');
end


end
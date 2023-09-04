function rmap = rasterMaps(x, y, signal, varargin)
%% rasterMaps
% Generate rastermap(s) of size [max(x_pos), max(y_pos)] of any type of 
% timeseries signal(s). 
%
% INPUT
%   --- Required
%   x (1 x N num array): The x_pos array can contain multiple
%       repetitions of the same x value, for instance [1 1 2 2 2 3 4 4]. In
%       this case, the mean of each indices in the signal input will be
%       averaged.
%   y (1 x N num array): The y_pos array will typically be a lap or
%       trial type of number, such that y_pos is something like [1 1 1 1 .... 2
%       2 2 2] etc.  
%   signal (M x N num array): This can be a matrix or array. If matrix is
%       given as input, the output rmap contains a rastermap for each row in the
%       matrix.
%
%   --- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters".
%       
% OUTPUT
%   rmap (Y x X x Z matrix): rmap will contain a rastermap for each row
%       in the signal input. If signal is a 1 x N array, rmap has only 2 dims 
%       (Y x X).
%
% EXAMPLES
%   This single function can be used to create a bunch of plots, for example:
%   rastermaps of running speed, image alignment, place cells, trial-by-trial
%   raster maps, licking profile to name a few.
%
%   1) PLACE CELL RASTERMAP
%       To create a rastermap of 100 roi signals at once, where we exclude
%       indices where the running speed is below 2 cm/s. 
%       The input we need is:
%       position: An array of length 1 to N, where each value is the position 
%           of the animal at that given sample of the input signal.
%       lap: An array of length 1 to N, where each value is the lap number
%           of the animal at that given sample of the input signal.
%       roiSignals: A 100 x N matrix, where each row in the matrix is the
%           signal of a roi.
%       runSpeed: A 1 x N array of the running speed at each sample of the
%           input signal.
%       
%       placeCellMaps = sb.generate.rasterMaps(position,lap,roiSignals,"exclude_inds",[runSpeed<2]); 
%
%   2) RUN SPEED RASTERMAP
%       To create a rastermap of the running speed: The minimum required input is
%       x_pos (which is a 1 x N array), y_pos (as a 1 x N array) and signal 
%       (which is a 1 x N array) containing the running speed. We also
%       tweak the smoothing parameter to be 10 inds. 
%
%       runrasterMap = sb.generate.rasterMaps(position,lap,runSpeed,"smoothing",10);
%
% Written by Andreas S Lande 2019

%% Default parameters
% Define default internal parameters used in the function
params = struct(); % The main container for all parameters
params.smoothing = 1; % Gaussian smoothing of the final rastermap
params.position_normalize = true; % Apply position normalization to each rastermap point
params.exclude_inds = zeros(1,length(x)); % Bool array where elements of true will result in the sample being excluded when generating the rastermap. This can be useful if for instance you want to create a rastermap of place cells but set a running threshold, or you want to exclude certain laps or trials

% Update parameters based on input
params = sb.helper.updateParameterStruct(params,varargin);

%% Check that inputs are correct
% Round X and Y pos inputs and make sure they are positive values
x = round(x);
y = round(y);
if min(x) < 0; x = x + abs(min(x)) + 1; elseif min(x) < 1; x = x + 1; end
if min(y) < 0; y = y + abs(min(y)) + 1; elseif min(y) < 1; y = y + 1; end

% Check that x_pos, y_pos and signal has the same length.
if ~isequal(size(x,2),size(y,2),size(signal,2),size(params.exclude_inds,2))
    warning("The inputs do not have the same length! Using smallest as reference.");
    minimum_length = min([size(x,2),size(y,2),size(signal,2),size(params.exclude_inds,2)]);
    x = x(:,1:minimum_length);
    y = y(:,1:minimum_length);
    signal = signal(:,1:minimum_length);
    params.exclude_inds = params.exclude_inds(:,1:minimum_length);
end

%% Create rastermap(s)
% Find the dimensions of the rastermap
x_max = max(x);
y_max = max(y);
n_maps = size(signal,1); % Number of signals to be used

% Specify number of samples
n_samples = length(x);

% Create rastermap
rmap = zeros(y_max, x_max, n_maps); % init rmaps to be created
posmap = zeros(y_max, x_max, n_maps); % init position maps to be used for normalize the rmaps to number of visits in each x, y position.

for n = 1:n_maps % For each map  
    for t = 1:n_samples % For each data sample
        if (params.exclude_inds(t) == 0) % Don't include excluded samples
            x_val = x(t);
            y_val = y(t);
            rmap(y_val,x_val,n) = rmap(y_val,x_val,n) + signal(n,t); % Add the signal to rmap
            posmap(y_val,x_val,n) = posmap(y_val,x_val,n) + 1; % Increment this X,Y position with 1 "visit".
        end
    end
end

% Create a position normalized rastermap
if params.position_normalize
    rmap = (rmap./posmap); 
end

% Smooth the rastermaps
if params.smoothing > 0
    rmap = smoothdata(rmap,2,'gaussian',params.smoothing);
end


end
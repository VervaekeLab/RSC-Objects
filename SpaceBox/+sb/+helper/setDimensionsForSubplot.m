function [x_dim,y_dim] =  setDimensionsForSubplot(n_axes)
%% setDimensionsForSubplot
% Find optimal dimensions when subplotting many images.
%
% INPUT
%   -- Required
%   n_axes: Number of elements to be plotted.
%
% OUTPUT
%   x_dim: Number of plots in x dimension.
%   y_dim: Number of plots in y dimension.
%
% Written by Andreas S Lande 2019

%% Find dimensions
if n_axes <= 5
    x_dim = n_axes;
    y_dim = 1;
elseif n_axes == 6
    x_dim = 3;
    y_dim = 2;
elseif n_axes == 8
    x_dim = 4;
    y_dim = 2;
elseif n_axes <= 10
    x_dim = 5;
    y_dim = 2;
elseif n_axes == 12
    x_dim = 4;
    y_dim = 3;
elseif n_axes <= 15
    x_dim = 5;
    y_dim = 3;
elseif n_axes <= 40
    x_dim = 10;
    y_dim = ceil(n_axes/10);
elseif n_axes <= 120
    x_dim = 15;
    y_dim = ceil(n_axes/15);
end



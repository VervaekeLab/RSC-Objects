line_width = 2; 

% Get the rastermap of running speed for each session
for s = 1:length(mData)    
    
    % Extract timeseries data
    pos = mData(s).sData.behavior.wheelPosDsBinned;
    lap = mData(s).sData.behavior.wheelLapDsBinned;
    run = mData(s).sData.behavior.runSpeedDs;
    lick = mData(s).sData.daqdata.lickSignal(mData(s).sData.daqdata.frameIndex);
    
    % Create rastermap of the running speed binned along the track
    running_rmap = sb.generate.rasterMaps(pos,lap,run,'smoothing',1);
    lick_rmap = sb.generate.rasterMaps(pos,lap,lick,'smoothing',1);
    
    % Find indices for laps when the landmark is on or off.
    on_laps = find([mData(s).sData.landmarks.motor.on(:,1) == 1]);
    
    if s == 1
        running_rmap_landmark_on = running_rmap(on_laps,:);        
        lick_rmap_landmark_on = lick_rmap(on_laps,:);
        
    else
        running_rmap_landmark_on = cat(1,running_rmap_landmark_on,running_rmap(on_laps,:));        
        lick_rmap_landmark_on = cat(1,lick_rmap_landmark_on,lick_rmap(on_laps,:));
    end
    
   
end

%% Plot individual runs and the average
figure(1); clf;
color_map = [0,0,0];
line_width = 1;
subplot(4,1,1:3)
% Plot speed when landmark is present
ops = struct();
ops.error = 'std';
ops.handle     = figure(1);
ops.alpha      = 0.1;
ops.line_width = line_width;
ops.color_area = color_map(1,:);
ops.color_line = color_map(1,:);
plot_areaerrorbar(running_rmap_landmark_on,ops);

hold on;

% Draw lines indicating the landmarks
line([32,32],[-100,100],'Color','r')
line([72,72],[-100,100],'Color','r')
line([92,92],[-100,100],'Color','k')

xticks([0,105])
xlim([0,105])
xticklabels({0,157})
title('Running speed');
ylabel('Running speed (cm/s)');
xlabel('Position (cm)');
ylim([0,58])
yticks([0,50])
set(gca,'FontSize',16);
box('off')

% Lick signal
lick_probability_on = nansum(lick_rmap_landmark_on,1);
lick_probability_on = lick_probability_on/sum(lick_probability_on);
lick_signal = nansum(lick_probability_on,1);
%lick_signal = normalize(lick_signal,'range',[0,20]);
subplot(4,1,4);
area(lick_signal,'FaceColor',[0,0,1],'FaceAlpha',0.4,'EdgeAlpha',0,'EdgeColor',[0,0,1],'LineWidth',line_width,'BaseValue',0);
hold on

xticks([0,105])
xlim([0,105])
xticklabels({0,157})
ylabel('Probability');
xlabel('Position (cm)');
%ylim([0,23])
%yticks([0,20])
set(gca,'FontSize',16);
box('off')


% Set figure size
set(gcf,'renderer', 'painters', 'Position', [100 200 800 450])


% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function plot_areaerrorbar(data, options)
    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'std';
    end
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = nanmean(data,1);
    data_std  = nanstd(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case '2std', error = data_std*2;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    end
    
    
    % Plotting the result
    figure(options.handle);
    x_vector = [options.x_axis', fliplr(options.x_axis')];
    patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
    'LineWidth', options.line_width);
    hold off;
    
    
    
end

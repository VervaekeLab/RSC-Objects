function axs = rasterMaps(rmaps, varargin)
%% rasterMaps 
% A quick plot of rastermaps, for example generated using 
% sb.generate.rasterMaps().
%
% INPUT
%   -- Required
%   rmaps: Y x X x N matrix of rastermap(s), where N is n number of rastermaps.
%       OR rmaps can also be a cell array where each element is a rmaps matrix
%           itself, like {rmaps1, rmaps2, rmaps3}. In this case each map will
%           be concatenated along the y dimension and plotted with separation
%           lines in between.
%       OR rmaps can be a single number n, in that case, a figure is created
%           with n blank rastermaps. When axs is taken as output, the specific
%           axes can therefor later be used to plot specific rastermaps in each
%           raster axes.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   axs: A axes array containing all axes in the plot.
%
% EXAMPLE
%   1) PLOT RASTERMAP OF DECONVOLVED SIGNAL PLACE CELLS
%       Because the function is automatically clipping minimum activity
%       to 1 by default, we set normalize_rastermap to true. Alternatively,
%       the parameter "color_map_minimum_max" can be lowered.
%
%       roi_rmaps = sb.generate.rasterMaps(position,laps,deconv_signal,'smoothing',6);
%       sb.plot.rasterMaps(roi_rmaps,'normalize_rastermap',true);
%
% Written by Andreas S Lande 2019


%% Default parameters
params = struct();
params.cmap = "jet"; % Colormap used for rastermaps
params.plot_average = true; % Enable plotting of average signal at bottom of rastermap
params.rmap_title = []; %num2cell(1:size(rmaps,3)); % Index number of each map
params.color_map_minimum_max = 0.7; % All maps have max value at minimum this value
params.color_map_max = []; % Max value for all colormaps. If empty, this is calculated for each raster map
params.show_max_number = true; % Show the max activity number for the rastermap in upper left corner
params.normalize_rastermap = false; % Normalize the whole rastermap to itself, range (0,1).
params.smoothing = 0; % Smoothing factor to be applied to rastermap. If set to 0 smoothing is omitted.
params.smoothing_average_plot = 0; % Number of bins to be used for smoothing the average plot if used. If set to 0 smoothing is omitted.
params.smoothing_method = "Gaussian"; % See MATLAB function smoothdata() for methods.
params.xlabel = ''; % Label in x dimension for each rastermap
params.ylabel = ''; % Label in y dimension for each rastermap
params.font_size = 12; % Size of text in plot
params.add_extra_plot = []; % An extra plot can be added, for example running speed profile. extra_plot must be of size 1 x size(rmaps,2).
params.add_lines_vertical = []; % x inds for adding vertical lines through plot
params.add_lines_horizontal = []; % y inds for adding horizontal lines through plot
params.colorbar = "off"; % on or off colorbar
params.axes = []; % Axes for each rastermap to be plotted.
params.add_on_color = 'w'; % The color choosen for the text and average plots added to the rastermaps
params.text_position_x = 5; % Pixel value for where the text is added
params.fig_size_x = 1000; % Size of figure in pixels
params.fig_size_y = 800; % Size of figure in pixels
params.correct_x_values_for_bin_size = 1.5; % Changes the values of x axis to fit the given bin size of the data.
params.cmap_max_type = 'peak'; % Specify how the max value of the cmap in imagesc should be specified. "peak" or "z_score"
params.max_plots_per_figure = 40; % Number of plots per figure created
params.gridAxesGapDims = [0.05,0.05]; % Gap size between each subplot in figure, x and y dim
params.markers = [];
params.markers_colors = [];
params.average_plot_max_method = 'compare'; % This sets the way the average plot max value is set. 'unique' will produce a unique average for each rastermap, while 'compare' will compare all maps towards each other.
params.super = true; % if super is set to true, a lot of default plotting is used. See below.

% Add super plots
if params.super
   params.smoothing = 3;
   params.add_lines_vertical = [32,72,92];
end

% Print all variables if "help" is only input 
try
    if rmaps == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

% Set params based on data
params = sb.helper.setParametersBasedOnRasterMapData(rmaps,params);

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);


%% If axes is given as input, clean current data in axes
if ~isempty(params.axes)
    for x = 1:length(params.axes)
        cla(params.axes(x));
    end
end

%% Find out if rmaps is cell array of rmaps

% If the input is a single number, generate a pseudo figure
if sum(size(rmaps)) == 2
    rmaps = zeros(10,10,rmaps);
    n_rows = size(rmaps,1);
    n_rows_each_map = n_rows;
    
else
    
    % If rmaps is a cell array; each element is treated as a unique rmaps that
    % shall be plotted on top of each other. The number of rows are thereby the
    % cumulative n rows for each of these maps.

    if iscell(rmaps) % Multiple maps to be merged
        n_rows = 0;
        n_rows_each_map = [];

        % For each cell array get number of rows
        for x = 1:length(rmaps)
           unique_map = rmaps{x};
           n_rows = n_rows + size(unique_map,1);
           n_rows_each_map(x) = size(unique_map,1);
        end

        % Merge all rmaps in cell array into single rmaps matrix
        unique_map = rmaps{1};
        newmaps = unique_map;
        for x = 2:length(rmaps)
            unique_map = rmaps{x};
            newmaps = cat(1,newmaps,unique_map);
        end

        % Update rmaps to be the merged newmaps matrix.
        rmaps = newmaps;

        % Set horizontal lines to each row of the map
        params.add_lines_horizontal = cumsum(n_rows_each_map);

    else % Only a single map exist
        n_rows = size(rmaps,1);
        n_rows_each_map = n_rows;
    end

end

%% Plot
n_maps = size(rmaps,3);
if n_maps < params.max_plots_per_figure
    params.max_plots_per_figure = n_maps;
end

number_of_figures = ceil(n_maps/params.max_plots_per_figure);
maps_left = 1:n_maps;

% Set the number of subfields to plot based on the number of maps
[x_dim,y_dim] =  sb.helper.setDimensionsForSubplot(params.max_plots_per_figure);
subplot_fields = [y_dim,x_dim];

% Create each figure needed to plot all ROIs
for figure_num = 1:number_of_figures
    
    % Find which rois to plot next
    if length(maps_left) > params.max_plots_per_figure
        maps_to_plot = maps_left(1:params.max_plots_per_figure);
        maps_left = maps_left(params.max_plots_per_figure+1:end);
    else
        maps_to_plot = maps_left(1:end);
        maps_left = [];
    end
    
    % Make figure
    if isempty(params.axes)
        figure('Renderer', 'painters', 'Position', [100 100 params.fig_size_x params.fig_size_y]);
        axs = multiPanelFigure.createAxesGrid(subplot_fields(1),subplot_fields(2),params.gridAxesGapDims);
    else
        axs = params.axes;
    end
    
    counter = 1;
    
    for map_index = maps_to_plot
        
        % Get current rastermap
        current_map = rmaps(:,:,map_index);
        
        % Normalize rastermap to itself
        if params.normalize_rastermap
            current_map = normalize(current_map,'range');
        end
        
        % Smooth the rastermap
        if params.smoothing
            current_map = smoothdata(current_map,2,params.smoothing_method,params.smoothing);
        end
        
        % Apply max min colormapping
        switch params.cmap_max_type
            case 'peak'
                cMin = 0; cMax = params.color_map_minimum_max;
                if max(max(current_map)) > params.color_map_minimum_max
                    cMax = max(max(current_map));
                end
                
            case 'z_score'
                z_val = 3;
                cMin = 0; cMax = params.color_map_minimum_max;
                if mean(current_map(:)) + z_val*std(current_map(:)) > params.color_map_minimum_max
                    cMax = mean(current_map(:)) + z_val*std(current_map(:));
                end
                
        end
           
        % Set max value to cmap if specified
        if ~isempty(params.color_map_max)
            cMax = params.color_map_max;
        end

        % Generate plot
        %axes(axs(counter));
        %subplot(subplot_fields(1),subplot_fields(2),counter);
        him = imagesc(current_map,[cMin,cMax]);
        him.Parent = axs(counter);
        axs(counter).CLim = [cMin,cMax];
        
        % Plot the max activity used for cMax
        if params.show_max_number
            text(axs(counter),params.text_position_x,round(n_rows/15),sprintf('%.2f',cMax),'Color',params.add_on_color,'FontSize',16,'FontWeight','bold');
        end
        
        % Set plot parameters
        %colormap(params.cmap);
        hold(axs(counter),'on');
        
        % Plot average
        if params.plot_average
            
            
            % Individual map laps
            map_start_inds = [1,cumsum(n_rows_each_map)+1];
            
            % Generate average for each map
            for x = 1:length(n_rows_each_map)
                average_map = [];
                average_map = nanmean(current_map(map_start_inds(x):map_start_inds(x+1)-1,:));
                
                average_map_max_value = nanmax(average_map);
                average_map = normalize(average_map,'range',[0,1]);
                average_map = average_map*(round(nanmean(n_rows_each_map))*0.25);
                
                if strcmp(params.average_plot_max_method,'compare')
                    average_map = average_map * (average_map_max_value / nanmax(nanmean(current_map)));
                end  
                
                % Smooth average plot
                if params.smoothing_average_plot
                    average_map = smoothdata(average_map,params.smoothing_method,params.smoothing_average_plot);
                end
                
                % Plot average at bottom of the rmap
                plot(axs(counter),-average_map+map_start_inds(x+1),'Color',params.add_on_color,'LineWidth',2)
            end
            
            
            
            
        end
        
        % Add vertical lines in plots
        if length(params.add_lines_vertical)>0
            % For each line index
            for x = 1:length(params.add_lines_vertical)
                line(axs(counter),[params.add_lines_vertical(x),params.add_lines_vertical(x)],[1,n_rows],'Color',[0.9,0.9,0.9],'LineWidth',2) 
            end
        end
        
        % Add horizontal lines in plots
        if length(params.add_lines_horizontal)>0
            % For each line index
            for x = 1:length(params.add_lines_horizontal)
                line(axs(counter),[1,size(rmaps,2)],[params.add_lines_horizontal(x),params.add_lines_horizontal(x)],'Color',params.add_on_color,'LineWidth',2) 
            end
        end
        
        % Add extra plot in rastermaps
        if length(params.add_extra_plot) > 0
    
            % Normalize extra plot to 20 % of rmaps size
            extra_plot = normalize(params.add_extra_plot,'range');
            extra_plot = extra_plot*(n_rows*0.2);
            
            % Plot average at bottom of the rmap
            plot(-extra_plot+n_rows,'Color',[0.7,0.7,0.7],'LineWidth',1.5)

        end
        
        % Add markers map
        % This is a star * map of signals over the map, for example to
        % indicate events somewhere in the rastermap
        if ~isempty(params.markers)
           
            for x = 1:length(params.markers)
               mm = params.markers{x};
               hold on; plot(axs(counter),mm(2,:),mm(1,:),'LineStyle','none','Marker','.','MarkerSize',8,'Color',params.markers_colors(x,:));
            end       
        end
          
        set(axs,'YDir','reverse');
        hold(axs(counter),'off');
        
        % Add colorbar
        if strcmp(params.colorbar,'on')
            colorbar;
        end
        
        % Set title to subplot
        if ~isempty(params.rmap_title)
            title(axs(counter),params.rmap_title{map_index},'FontSize',17);
        else
            title(axs(counter),map_index,'FontSize',17);
        end
        
        % Set plot texts
        xlabel(params.xlabel);
        ylabel(params.ylabel);
        set(axs(counter),'FontSize',params.font_size);
        
        % Set other plot parameters
        set(axs(counter),'YLim',[1,size(current_map,1)])
       
        xtix = get(axs(counter), 'XTick');
        set(axs(counter), 'XTick',xtix, 'XTickLabel',xtix*params.correct_x_values_for_bin_size);
        set(axs(counter),'FontSize',18);
        
        counter = counter + 1;
    end

    % Apply axes properties to all   
    %set(axs,'XTick',[]);
    %set(axs,'YTick',[]);
    colormap(params.cmap);
    
end

end

%% Create example of the task structure and show example cell responses

% Parameters
landmark_line_width = 1;
session_to_plot = 3;
signal_type = 'dffSubtractedNpil';

% ---- Plot the task structure
for q = 1
    
    figure(1); clf;
    
    % Show the structure of the passive flyby task
    landmark_color = [239,35,60]/255;
    water_color = [42,150,255]/255;
    lick_color = [0.5,0.5,0.5];
    figure(1); clf;
    
    % ---- Enganged
    subplot(3,6,1:3);
    start_ind = 1075;
    end_ind = 1890;
    window_size = end_ind-start_ind;
    
    % Motor signal
    motor_signal = normalize(mData(session_to_plot).sData.daqdata.left_motor_speed(mData(session_to_plot).sData.daqdata.frameIndex),'range',[-10,1000]);
    
    % Draw a simple line instead of the actual motor signal
    motor_signal = motor_signal(start_ind:end_ind);
    motor_signal(motor_signal>1) = 1;
    motor_signal(motor_signal<0) = 0;
    motor_signal = [0,diff(motor_signal)];
    motor_signal(motor_signal<0) = 0;
    motor_signal = find(motor_signal);
    
    % Draw line at each motor onset
    for i = 1:length(motor_signal)
        line([motor_signal(i),motor_signal(i)],[-1000,1000],'Color','r','LineWidth',landmark_line_width)
    end

    hold on
    
    % Lick signal
    lick_signal = normalize(mData(session_to_plot).sData.daqdata.lickSignal(mData(session_to_plot).sData.daqdata.frameIndex),'range',[-10,1000]);
    
    % Draw tick lick signal
    lick_signal = lick_signal(start_ind:end_ind);
    lick_signal(lick_signal<1) = nan;
    lick_signal(lick_signal>1) = 50;
    plot(lick_signal,'Marker','.', 'MarkerSize',18,'Color',[0.5,0.5,0.5])
    
    % Water signal
    water_signal = normalize(mData(session_to_plot).sData.daqdata.waterValveSignal(mData(session_to_plot).sData.daqdata.frameIndex),'range',[-10,1000]);
    
    % Draw simplified water signal
    water_signal = water_signal(start_ind:end_ind);
    water_signal(water_signal>1) = 1;
    water_signal(water_signal<0) = 0;
    water_signal = [0,diff(water_signal)];
    water_signal(water_signal<0) = 0;
    water_signal = find(water_signal);
    
    % Draw line at each motor onset
    for i = 1:length(water_signal)
        line([water_signal(i),water_signal(i)],[-1000,1000],'Color',water_color,'LineWidth',landmark_line_width)
    end
    
    % Running speed
    plot(mData(session_to_plot).sData.behavior.runSpeedDs(start_ind:end_ind),'Color','k','LineWidth',1)
    
    % Properties
    ylim([-5,55])
    yticks([0,50])
    xlim([0,end_ind-start_ind])
    ylabel(sprintf('Running speed\n(cm/s)'));
    title('Enganged');
    set(gca,'FontSize',16);
    xticks([]);
     
    % ---- Passive (immobile)
    subplot(3,6,13:15);
    start_ind = 1000;
    end_ind = start_ind + window_size;
    
    plot(mData(session_to_plot+1).sData.behavior.runSpeedDs(start_ind:end_ind),'Color','k','LineWidth',1)
    hold on;
    
    % Motor signal
    motor_signal = normalize(mData(session_to_plot+1).sData.daqdata.left_motor_speed(mData(session_to_plot+1).sData.daqdata.frameIndex),'range',[-10,1000]);
    
    % Draw a simple line instead of the actual motor signal
    motor_signal = motor_signal(start_ind:end_ind);
    motor_signal(motor_signal>1) = 1;
    motor_signal(motor_signal<0) = 0;
    motor_signal = [0,diff(motor_signal)];
    motor_signal(motor_signal<0) = 0;
    motor_signal = find(motor_signal);
    
    % Draw line at each motor onset
    for i = 1:length(motor_signal)
        line([motor_signal(i),motor_signal(i)],[-1000,1000],'Color',landmark_color,'LineWidth',landmark_line_width)
    end
    
    ylim([-5,55])
    yticks([0,50])
    xlim([0,window_size])
    ylabel(sprintf('Running speed\n(cm/s)'));    
    title('Passive stimulus (immobile)');
    xlabel('Time (seconds)');
    xticks([1,window_size])
    xticklabels({0,26})
    
    % Draw a line of length 3 seconds to use as time indicator
    line([700,793],[40,40],'LineWidth',1,'Color','k')

    set(gca,'FontSize',16);
    
    
    % ---- Passive (running)
    subplot(3,6,7:9);
    start_ind = 18000;
    end_ind = start_ind + window_size;
    
    % Motor signal
    motor_signal = normalize(mData(session_to_plot+1).sData.daqdata.left_motor_speed(mData(session_to_plot+1).sData.daqdata.frameIndex),'range',[-10,1000]);
    
    % Draw a simple line instead of the actual motor signal
    motor_signal = motor_signal(start_ind:end_ind);
    motor_signal(motor_signal>1) = 1;
    motor_signal(motor_signal<0) = 0;
    motor_signal = [0,diff(motor_signal)];
    motor_signal(motor_signal<0) = 0;
    motor_signal = find(motor_signal);
    
    % Draw line at each motor onset
    for i = 1:length(motor_signal)
        line([motor_signal(i),motor_signal(i)],[-1000,1000],'Color',landmark_color,'LineWidth',landmark_line_width)
    end
    
    hold on;

    % Running speed
    plot(mData(session_to_plot+1).sData.behavior.runSpeedDs(start_ind:end_ind),'Color','k','LineWidth',1)
   
    ylim([-5,55])
    xlim([0,window_size])
    yticks([0,50])
    ylabel(sprintf('Running speed\n(cm/s)'));
    title('Passive stimulus (running)');
    xticks([]);
    set(gca,'FontSize',16);
    
      
    % Set correct size of the figure
    set(gcf,'renderer', 'painters', 'Position', [-1200,300,1000,500])
end

% ---- Plot example cell responses
% Shaded error plot options
options.handle     = figure(1);
options.color_area = [150 150 150]./255;    % Blue theme
options.color_line = [ 0 0 0]./255;
options.alpha      = 0.8;
options.line_width = 0.5;
options.error      = 'sem';

% Example Cell #1 - Landmark cell that is NOT respondin to passive stimulus
for q = 2
   
    % Landmark cell that is NOT respondin to passive stimulus
    c_analysis = 31 %38; % This is the cell index in the analysis struct.
    
    % Actual session and roi index
    s = analysis(c_analysis).session
    c = analysis(c_analysis).ind;
    
    % Find the dff signal following the touch of the motor in engaged task
    dff = mData(s).sData.imdata.roiSignals(2).(signal_type)(c,:);

    % Find each bin where the mouse enters position 50 cm.
    pos = mData(s).sData.behavior.wheelPosDs;
    
    % Find each bin where the animal passes the landmark
    
    % Set some init vars
    bins_at_landmark = [];
    landmark_passed = true;
    landmark_position = 48; % the estimated position of the first position where the animal is in contact with the landmark ( in cm. )
    i = 0;
    
    while i < length(pos)
        i = i+1;
        if pos(i) > landmark_position
            
            if ~landmark_passed
                bins_at_landmark = [bins_at_landmark,i];
                landmark_passed = true;
                i = i+20; % Do a jump so that we don't get multiple values in case there are more interactions with the landmark
            end
        else
            landmark_passed = false;
        end
        
    end
    
    % Get the dff at each landmark passing
    window_size = 62; % time bins before and after the landmark interaction
    window_n_bins = length([1-window_size:1+window_size]);
    
    dff_at_landmark = nan(length(bins_at_landmark)-50,window_n_bins);
    for i = 1:length(bins_at_landmark)-50
        try
            dff_at_landmark(i,:) = dff(bins_at_landmark(i)-window_size:bins_at_landmark(i)+window_size);
        end
    end

    % Add the plot to the existing figure 1
    % Enhanced
    figure(1);
    subplot(3,6,4);
    plot_areaerrorbar(dff_at_landmark(2:end,:),options)
    line([window_size,window_size],[-1,3],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F')
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    title('Cell #1')
    set(gca,'FontSize',16)

    % Passive (immobile)
    figure(1);
    subplot(3,6,16);
    plot_areaerrorbar(analysis(c_analysis).still.rmap(:,62:186),options)
    line([window_size,window_size],[-1,3],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F')
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
        xlabel('Time (seconds)');
    set(gca,'FontSize',16)

    % Passive (running)
    figure(1);
    subplot(3,6,10);
    plot_areaerrorbar(analysis(c_analysis).moving.rmap(:,62:186),options)
    line([window_size,window_size],[-1,3],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F')
    xlim([0,window_n_bins])
    
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    set(gca,'FontSize',16)
    
    
end

% Example Cell #2 - Cell that responds to passive stimulus
for q = 3
   
    % Landmark cell that responds to passive stimulus
    c_analysis = 14; % This is the cell index in the analysis struct.
    
    % Actual session and roi index
    s = analysis(c_analysis).session
    c = analysis(c_analysis).ind;
    
    % Find the dff signal following the touch of the motor in engaged task
    dff = mData(s).sData.imdata.roiSignals(2).(signal_type)(c,:);

    % Find each bin where the mouse enters position 50 cm.
    pos = mData(s).sData.behavior.wheelPosDs;
    
    % Find each bin where the animal passes the landmark
    
    % Set some init vars
    bins_at_landmark = [];
    landmark_passed = true;
    landmark_position = 48; % the estimated position of the first position where the animal is in contact with the landmark ( in cm. )
    i = 0;
    
    while i < length(pos)
        i = i+1;
        if pos(i) > landmark_position
            
            if ~landmark_passed
                bins_at_landmark = [bins_at_landmark,i];
                landmark_passed = true;
                i = i+20; % Do a jump so that we don't get multiple values in case there are more interactions with the landmark
            end
        else
            landmark_passed = false;
        end
        
    end
    
    % Get the dff at each landmark passing
    window_size = 62; % time bins before and after the landmark interaction
    window_n_bins = length([1-window_size:1+window_size]);
    
    dff_at_landmark = nan(length(bins_at_landmark)-50,window_n_bins);
    for i = 1:length(bins_at_landmark)-50
        try
            dff_at_landmark(i,:) = dff(bins_at_landmark(i)-window_size:bins_at_landmark(i)+window_size);
        end
    end

    % Add the plot to the existing figure 1
    % Enhanced
    figure(1);
    subplot(3,6,5);
    plot_areaerrorbar(dff_at_landmark(2:end,:),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    title('Cell #2')
    set(gca,'FontSize',16)
    
    % Passive (immobile)
    figure(1);
    subplot(3,6,17);
    plot_areaerrorbar(analysis(c_analysis).still.rmap(:,62:186),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    xlabel('Time (seconds)');

    set(gca,'FontSize',16)

    
    % Passive (running)
    figure(1);
    subplot(3,6,11);
    plot_areaerrorbar(analysis(c_analysis).moving.rmap(:,62:186),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');

    
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    set(gca,'FontSize',16)
    
    
    
    
end

% Example Cell #3 - Cell that responds only to passive stimulus
for q = 4
   
    % Landmark cell that responds to passive stimulus
    c_analysis = 11; % This is the cell index in the analysis struct.
    
    % Actual session and roi index
    s = analysis(c_analysis).session
    c = analysis(c_analysis).ind;
    
    % Find the dff signal following the touch of the motor in engaged task
    dff = mData(s).sData.imdata.roiSignals(2).(signal_type)(c,:);

    % Find each bin where the mouse enters position 50 cm.
    pos = mData(s).sData.behavior.wheelPosDs;
    
    % Find each bin where the animal passes the landmark
    
    % Set some init vars
    bins_at_landmark = [];
    landmark_passed = true;
    landmark_position = 48; % the estimated position of the first position where the animal is in contact with the landmark ( in cm. )
    i = 0;
    
    while i < length(pos)
        i = i+1;
        if pos(i) > landmark_position
            
            if ~landmark_passed
                bins_at_landmark = [bins_at_landmark,i];
                landmark_passed = true;
                i = i+20; % Do a jump so that we don't get multiple values in case there are more interactions with the landmark
            end
        else
            landmark_passed = false;
        end
        
    end
    
    % Get the dff at each landmark passing
    window_size = 62; % time bins before and after the landmark interaction
    window_n_bins = length([1-window_size:1+window_size]);
    
    dff_at_landmark = nan(length(bins_at_landmark)-50,window_n_bins);
    for i = 1:length(bins_at_landmark)-50
        try
            dff_at_landmark(i,:) = dff(bins_at_landmark(i)-window_size:bins_at_landmark(i)+window_size);
        end
    end

    % Add the plot to the existing figure 1
    % Enhanced
    figure(1);
    subplot(3,6,6);
    plot_areaerrorbar(dff_at_landmark(2:end,:),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    title('Cell #3')
    set(gca,'FontSize',16)
    
    % Passive (immobile)
    figure(1);
    subplot(3,6,18);
    plot_areaerrorbar(analysis(c_analysis).still.rmap(:,62:186),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    xlabel('Time (seconds)');

    set(gca,'FontSize',16)

    
    % Passive (running)
    figure(1);
    subplot(3,6,12);
    plot_areaerrorbar(analysis(c_analysis).moving.rmap(:,62:186),options)
    line([window_size,window_size],[-1,1],'Color',landmark_color,'LineWidth',landmark_line_width)
    ylim([-0.1,0.5])
    yticks([0,0.5])
    ylabel('dF/F');

    
    xlim([0,window_n_bins])
    xticks([0, 62, 124])
    xticklabels({-2,0,2})
    set(gca,'FontSize',16)
    
    
    
    
end




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

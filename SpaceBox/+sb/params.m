function p = params(name,set_value)
%% PARAMS
% This is used to set parameters for SpaceBox functions to extract it. The 
% INPUT:
%   name: name of function
%   set_value: (optional), if set_value is given, the current value of the
%   parameter will be overwritten in parameter file.
%
% OUTPUT:
%   p: output value of the parameter asked for
%
% WB Andreas S Lande | Vervaeke Lab 2019
% 
% 


% CUE PARAMETERS (cm)
% The cue start points are set to 1.5 cm in front of mouse nose. 
sandpaper_start = [25.5,115.5]; % Nose is at 28.5 and 119.5.
sandpaper_end = [33,123]; 
velcro_start = 6; % OLD MOUSE: 4.5; % Nose is at 7 cm.
velcro_end = 10; % old mouse %10.5;
glue_spikes_start = [58,119];  % old mouse [48,84]; % Nose is at 51 cm and 87 cm
glue_spikes_end = [61,122];%[55.5,91.5];
        

% PARAMETERS TO SET
switch name 
    case 'bin_size_wheel' % The bin size in cm for each bin the wheel linear track is divided into. Make sure this corresponds with bins_n_wheel.
        p = 1.5;
        
    case 'bins_n_wheel' % Number of bins the wheel linear track is to be divided into. Make sure this corresponds with bin_size_wheel.
        p = 105;
        
    case 'cue_sandpaper' % Position of sandpaper cues
        % Sandpaper extends from 30 to 35 and from 121 to 126 cm position
        % When measured from mouse nose, they start at 28.5 and 119.5 cm
        
        sandpaper = zeros(1,105);
        sandpaper(round(sandpaper_start(1)/sb.params('bin_size_wheel')):round(sandpaper_end(1)/sb.params('bin_size_wheel'))) = 1;
        sandpaper(round(sandpaper_start(2)/sb.params('bin_size_wheel')):round(sandpaper_end(2)/sb.params('bin_size_wheel'))) = 1;
        p = sandpaper;
        
    case 'cue_sandpaper_cm' % Position of sandpaper cues in cm
        % Sandpaper extends from 30 to 35 and from 121 to 126 cm position
        % When measured from mouse nose, they start at 28.5 and 119.5 cm
        
        sandpaper = zeros(1,158);
        sandpaper(round(sandpaper_start(1)):round(sandpaper_end(1))) = 1;
        sandpaper(round(sandpaper_start(2)):round(sandpaper_end(2))) = 1;
        p = sandpaper;
        
        
    case 'cue_velcro' % Position of velcro cues
        % Velcro extends from 8.5 to 12.5 cm on wheel in relation to reward
        velcro = zeros(1,105);
        velcro(round(velcro_start/sb.params('bin_size_wheel')):round(velcro_end/sb.params('bin_size_wheel'))) = 1;
        p = velcro;
        
    case 'cue_velcro_cm' % Position of velcro cues
        % Velcro extends from 8.5 to 12.5 cm on wheel in relation to reward
        velcro = zeros(1,158);
        velcro_start = 6.5;
        velcro_end = 12.5;
        velcro(round(velcro_start):round(velcro_end)) = 1;
        p = velcro;
        
    case 'cue_glue_spikes' % Position of spike glue cues
        % Glue spikes extend from 52.5 to 56.5 and 88.5 to 92.5 cm.
        glue_spikes = zeros(1,105);
        glue_spikes(round(glue_spikes_start(1)/sb.params('bin_size_wheel')):round(glue_spikes_end(1)/sb.params('bin_size_wheel'))) = 1;
        glue_spikes(round(glue_spikes_start(2)/sb.params('bin_size_wheel')):round(glue_spikes_end(2)/sb.params('bin_size_wheel'))) = 1;
        p = glue_spikes;
        
    case 'cue_glue_spikes_cm' % Position of spike glue cues
        % Glue spikes extend from 52.5 to 56.5 and 88.5 to 92.5 cm.
        glue_spikes = zeros(1,158);
        glue_spikes(round(glue_spikes_start(1)):round(glue_spikes_end(1))) = 1;
        glue_spikes(round(glue_spikes_start(2)):round(glue_spikes_end(2))) = 1;
        p = glue_spikes;
        
    case 'font_size' % Font size of plots generated by SpaceBox
        p = 15;
        
    case 'fov_rotation' % Rotation of FOV used for showing the FOV image the correct way up.
        p = -90;
    case 'peak_specific_reliability' % This is either TRUE or FALSE. If true, a place field is only included in the peak count if it fullfills the place field reliability criteria. If this is zero, only the major field has to go through this test.
        p = 0; 
    case 'place_field_activity_ratio_criteria' % In-field vs out-of-field activity ratio criteria as specified for place cell detection by Mao et al (2017).
        p = 2;
    
    case 'place_field_max_peak_width' % Max place field width size in cm.
        p = 150;
    
    case 'place_field_min_peak_distance' % Minimum place field peak distance (between multiple peaks).
        p = 40;
    
    case 'place_field_min_peak_width' % Minimum place field width size in cm.
        p = 2;
           
    case 'place_field_reliability_criteria' % Number of laps the peak of the cells firing has to be within the averaged place field. Given as percentage / 100. So 0.5 is 50% reliability.
        p = 0.2;
        
    case 'reward_offset'
        p = 0; % bins
       
    case 'roi_signal_type' % Set the signal within roiSignals(ch). you wish to be used for plotting and analysis
        p = 'dffSubtractedNpil'; %'deconv';
        %p = 'deconv';
    case 'run_threshold'
        p = 2; % cm. This is the minimum speed required for a time sampled bin to be included in producing activity heat maps etc. Basically, each time bin where the animal is running at a less speed than this will be excluded.
        
    case 'smooth_deconvolve' % Smooth the deconvolved signal before plotting?
        p = 1;
        
    case 'smoothing' % Smooth the plots? Gaussian is used. Set to 0 if you dont want it. The number correspond to number of bins one the linear track.
        p = 3;
        
    case 'wheel_length' % Circumference of the running wheel.
        p = 157.5;
        

end
       

end


% % Load parameter file
% file_num = 1;
% if file_num == 1
%     file_path = fullfile(getPathToDir('parameters'),'parameters.xlsx');
% else
%     file_path = fullfile(getPathToDir('parameters'),'parameters2.xlsx');
% end
% 
% [~,~,rawdata] = xlsread(file_path);
% 
% 
% % Based on number of inputs
% if nargin > 1 % Set value
%     newdata = rawdata;
%     for x = 1:size(rawdata,1)
%         if strcmp(lower(rawdata{x,1}), lower(name))
%            newdata{x,2} = set_value;
%         end
%     end
%     
%     xlswrite(file_path,newdata);
%     
% else % Get value
%     
%     for x = 1:size(rawdata,1)
%         if strcmp(lower(rawdata{x,1}), lower(name))
%             p = rawdata{x,2};
%             % Convert to number if needed
%             try
%                 if strcmp(rawdata{x,3},'num')
%                    p = str2num(p); 
%                 end
%             end
%         end
%     end
%     
%     
% end
% 
% 
% 
% 
% 
% end

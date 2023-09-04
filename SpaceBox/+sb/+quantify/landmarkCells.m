function db = landmarkCells(db,varargin)
%% landmarkCells
% Quantifies the number of landmark cells in each recording.
%
% INPUT
%   -- Required
%   db: Database.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   db: Updated database.
%
% Written by Andreas S Lande 2019

% Params
params.raster_signal = 'deconv';
params.show_info = true; % Plot visualization of the quantification

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

%% Quantify
% For each session
for session = 1:length(db)
    
    % Init var
    landmark_cells = [];
    
    try
        % Get place cells for session
        place_cells = db(session).cells.place_cells_deconv;

        % Grab rastermaps for this session
        rmaps = db(session).roimaps.(params.raster_signal);

        % Grab landmark on laps
        landmarks_on = db(session).motor_on_laps;

        % Use only laps where the motor cues are on and Use only place cell rastermaps
        rmaps = rmaps(landmarks_on,:,place_cells);

        % Find landmark positions
        landmarks = round(unique([db(session).motor_landmarks.left,db(session).motor_landmarks.right])/1.5);
        % Add accepted offset surrounding all landmark positions on wheel
        accepted_offset_before = 5;
        accepted_offset_after = 8;
        lm = [];
        for l = 1:length(landmarks)
            lm{l} = landmarks(l)-accepted_offset_before:landmarks(l)+accepted_offset_after;
        end
        landmarks = lm;

        % For each roi in rastermap do analysis
        for r = 1:size(rmaps,3)

            % Estimate place cell statistics
            [~,stats] = sb.classify.placeTunedRoisSimple(rmaps(:,:,r),'min_peak_distance',10);

            in_field_inds = [stats.in_field_inds{:}];


            % Check that all landmarks are represented in the in_field_inds
            is_represented = 0;
            for l = 1:length(landmarks)
                if sum(ismember(landmarks{l},in_field_inds))
                   is_represented = is_represented + 1; 
                end
            end

            % If both landmarks are represented, add the cell number to landmark_cells
            if is_represented == 2
                landmark_cells = [landmark_cells,place_cells(r)];
            end  

        end
    end
    
    db(session).quantified.landmark_cells = landmark_cells;
    
end

% Create a plot showing the number of landmark cells across sessions
if params.show_info
    
    figure(7291299);
    clf;
    
    % Two subplots will be used
    subplot(1,4,1:3);
    
    % Grab percentage of landmark cells for each session
    percentage_landmark_cells = [];
    for s = 1:length(db)
        try
            percentage_landmark_cells(s) = length(db(s).quantified.landmark_cells) / length(db(s).cells.active)*100;
        catch
            percentage_landmark_cells(s) = nan;
        end
    end
    
    % Plot each session as a barplot
    bar(percentage_landmark_cells);
    hold on
    line([0,length(db)+1],[nanmean(percentage_landmark_cells),nanmean(percentage_landmark_cells)],'Color','r');
    ylabel('% of active cells');
    xlabel('Database ID');
    xticks(1:length(db));
    xticklabels([db(:).database_number]);
    title('Landmarks for each session');
    
    set(gca,'FontSize',15);
    
    % Plot the distribution
    subplot(1,4,4);
    histogram(percentage_landmark_cells,ceil(max(percentage_landmark_cells)));
    title('Landmark cell distribution');
    ylabel('# sessions');
    xlabel('% landmark cells')
    set(gca,'FontSize',15);
    

end


end
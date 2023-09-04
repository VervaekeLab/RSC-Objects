function [onsets_still, onsets_moving] = findLandmarkFlybyOnsets(run,motor,show_plot)
% Find the index for each sample of where a motor flyby occurs while the
% animal is sitting still before (and after)
%
% INPUT
%   run: Run signal
%   motor: Motor speed signal
%
%
% OUTPUT
%   onsets: Sample index for where a sample period starts.

if nargin < 3
    show_plot = false;
end

fs = 31; % Imaging frame rate
still_before = fs * 2; % fs * number of seconds the animal should be still before the motor onset
still_after = fs * 2; % fs * number of seconds the animal should be still after the motor onsets
run_threshold = 2; % cm/s

% Flatten the motor signal to be 1 where the motor starts moving and 0
% elsewhere
motor(motor>0) = 1;
motor_starts = diff(motor);
motor_starts(motor_starts<0) = 0;

% Find each sample index where the motor starts
motor_starts = find(motor_starts);

onsets_still = [];
onsets_moving = [];

% For each motor start search for periods where the animal sits still
for s = motor_starts
   try
        onset_window = run(s-still_before:s+still_after);

        % Add the sample index to the list if the animal does not move more
        % than allowed run threshold in the window.
        if max(abs(onset_window)) < run_threshold
            onsets_still = [onsets_still,s];
        else
            onsets_moving = [onsets_moving, s];
        end
   end
end


if show_plot
    % Plot all
    figure; clf;
    subplot(1,2,1);
    for x = 1:length(onsets_still)
        plot(run(onsets_still(x)-still_before:onsets_still(x)+still_before));
        hold on;
    end
    title(sprintf('Sitting still: %i',length(onsets_still)));
    y_lims = get(gca,'YLim');


    subplot(1,2,2);
    for x = 1:length(onsets_moving)
        plot(run(onsets_moving(x)-still_before:onsets_moving(x)+still_before));
        hold on;
    end

    title(sprintf('Moving: %i',length(onsets_moving)));

    % Set correct y lims
    y_lims = [y_lims,get(gca,'YLim')];
    subplot(1,2,1);
    set(gca,'YLim',[min(y_lims)-2,max(y_lims)+3]);
    subplot(1,2,2);
    set(gca,'YLim',[min(y_lims)-2,max(y_lims)+3]);
    end
end
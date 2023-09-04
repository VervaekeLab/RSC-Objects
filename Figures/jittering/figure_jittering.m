

dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/Nora/jittering/';

% Recordings to include
data.sessionIDs = { 'm1551-20230308-03', ... # Landmarks removed after 41 laps. Ran 96 laps total
                    'm1551-20230309-01', ... # Landmarks removed after 41 laps. Ran 136 laps total
                    'm1553-20230308-01', ... # Landmarks removed after 30 laps. Ran 44 laps total
                    'm1553-20230309-01', ...
                    'm1554-20230308-01',...
                    'm1554-20230309-01'
                    };
                

[mData,data] = rsc_jitter(dirct,data);

plot_sequence;

plot_bar;







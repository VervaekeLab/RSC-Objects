%% FIGURE C - When does the cell start to increase its firing preceding the landmark? 
% POSITION AXIS
predictive_cell_start_bin = ones(1,size(predictive_rmaps,3)-1)*105;
predictive_cell_max_bin = ones(1,size(predictive_rmaps,3)-1)*105;

% Threshold is 95th percentile of the averaged map
for c = 1:size(predictive_rmaps,3)
    avg_map = nanmean(predictive_rmaps(:,1:105,c));
    [~,max_bin] = nanmax(avg_map);
    predictive_cell_max_bin(c) = max_bin;
    threshold = prctile(avg_map,95);
    
    for b = 1:length(avg_map)
       if avg_map(1,b)>threshold
           predictive_cell_start_bin(c) = b;
           break;
       end
    end
        
end

predictive_cell_start_bin(predictive_cell_start_bin<33) = predictive_cell_start_bin(predictive_cell_start_bin<33) - 32;
predictive_cell_start_bin(predictive_cell_start_bin>33) = predictive_cell_start_bin(predictive_cell_start_bin>33) - 72;

prediction_start = sort(predictive_cell_start_bin*1.5);
figure(3); 
clf;
subplot(1,3,1);
histogram(prediction_start,'BinEdges',-50:5:0,'Normalization','probability');
ylim([0,0.25])
yticks([0.1,0.2,0.3])
yticklabels([10,20,30])
ylabel('% of cells');
xlabel('Distance (cm)')
set(gca,'FontSize',18)

% TIME AXIS
window_size = 31;
n_cells = length(predictive_cells_inds);

predictive_response = [];
    
for c = 1:n_cells
    
    % Find the indices when the mouse is entering into position 32 bin
    at_first_landmark = mData(predictive_cells_session(c)).sData.behavior.wheelPosDsBinned == 32;

    % Clean it up by setting bins where at_first_landmark has multiple 1's
    % after each other.
    at_first_landmark_clean = zeros(1,length(at_first_landmark));
    x = 1;
    while x < length(at_first_landmark)

        if at_first_landmark(x) == 1
            at_first_landmark_clean(x) = 1;
            x = x+5;
        else
            x = x+1;
        end

    end

    at_first_landmark_ind = find(at_first_landmark_clean);
    
    roi_response = mData(predictive_cells_session(c)).sData.imdata.roiSignals(2).deconv(predictive_cells_inds(c),:);
    
    for x = 1:length(at_first_landmark_ind)
        
        try
            predictive_response(x,:,c) = roi_response(at_first_landmark_ind(x)-window_size+1:at_first_landmark_ind(x));
        end
        
    end
    
end

for s = 1:length(mData)
    no_obj_cells(s) = length(unique([mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind]));
    prct_predictive_cells_session(s) = 100*length(find(predictive_cells_session == s))./no_obj_cells(s);
end

mean_prct_predictive_cells =  nanmean(prct_predictive_cells_session);
sem_prct_predictive_cells =  nanstd(prct_predictive_cells_session)/sqrt(length(mData));

print_prct_cells = 'There are mean +/- sem: %f +/- %f predictive cells across sessions \n';
fprintf(print_prct_cells ,mean_prct_predictive_cells,sem_prct_predictive_cells )

% Find the first ind where the there is a increase in the deconv signal
first_rise = [];

for c = 1:size(predictive_response,3)
    
    average_response = nanmean(predictive_response(:,:,c));
    threshold = prctile(average_response,95);
    
    first_rise(c) = 0;
    
    for t = 1:length(average_response)
        if average_response(t)>threshold
            if first_rise(c) == 0
                first_rise(c) = t;
            end
        end    
    end
    
end

figure(3); 
subplot(1,3,2)
histogram(first_rise,10,'Normalization','probability')
xlim([0,window_size])

yticks([0.1,0.2,0.3])
yticklabels([10,20,30])
ylabel('% of cells')
ylim([0,0.25])
xlabel('Time (seconds)')
xticks([1,window_size])
xticklabels({num2str(-window_size*(1/31)),'0'})
set(gca,'FontSize',18)



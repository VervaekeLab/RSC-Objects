%% Show FOV with all tactile cells and place cells
s = 7;
lcs = [mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind];
pcs = mData(s).rmaps.pcs.ind;

classes = zeros(1,size(mData(s).rmaps.dff,3));
classes(lcs) = 1;
classes(pcs) = 2;
blank_image = ones(size(mData(s).sData.imdata.meta.fovImageAvg,1),size(mData(s).sData.imdata.meta.fovImageAvg,2))*255;
blank_image(1,:) = 0;
blank_image(498,:) = 0;
blank_image(:,1) = 0;
blank_image(:,480) = 0;

sb.plot.roiColoredByTag(mData(s).sData.imdata.roiArray,classes,'fov_image',blank_image)


%% Plot the distribution of max amplitude for place vells vs landmark cells
% -------- Max amplitude
amp_pcs = [];
amp_lcs = [];
pctile = 99;
for s = 1:length(mData)
    
   
    pcs = mData(s).rmaps.pcs.ind;
    dff = mData(s).sData.imdata.roiSignals(2).dffSubtractedNpil(pcs,:);
    
    for r = 1:size(dff,1)
        amp_pcs = [amp_pcs,nanmax(dff(r,:))];
    end
    
    lcs = [mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind];
    dff = mData(s).sData.imdata.roiSignals(2).dffSubtractedNpil(lcs,:);
    
    for r = 1:size(dff,1)
        amp_lcs = [amp_lcs,nanmax(dff(r,:))];
    end
    
    
end

% Plot histograms normalized to probability
figure(1);clf;
subplot(1,2,1);
h = histogram(amp_pcs,'FaceColor','k','Normalization','probability','BinEdges',[0:0.2:4],'FaceAlpha',1,'FaceColor',[51,51,51]/255);
hold on;
h = histogram(amp_lcs,'FaceColor','r','Normalization','probability','BinEdges',[0:0.2:4],'FaceAlpha',0.9,'FaceColor',[239,35,60]/255);

title('dF/F max amplitude')
% Set labels
xlabel('dF/F');
ylabel('% of cells')
yticks([0,0.05,0.1,0.15])
yticklabels({0,5,10,15})
ylim([0,0.18])
xticks([0,4])

% Set limits and params
xlim([0,4])
set(gca,'FontSize',16);

% Plot the distribution of field width for place cells and landmark cells

% Landmark cells
lcs = [];
for s = 1:length(mData)
        
    if isempty(lcs)
        lcs = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.deconv_motor_on);
        if isempty(lcs)
            lcs = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on);
        else
            lcs = [lcs;sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on)];
        end
    else
        lcs = [lcs;sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.deconv_motor_on)];
        if ~isempty(mData(s).rmaps.tcs.ind)
            lcs = [lcs;sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on)];
        end
    end
    
end


lcs_field_width = [];

% Find field width for each cell at the first landmark position
for c = 1:size(lcs,1)
    response = lcs(c,10:65);
    [~,max_response_position(c)] = nanmax(response);
    threshold = nanmean(response);
    
%     figure(1); clf;
%     plot(response);
%     hold on;
%     line([max_response_position(c),max_response_position(c)],[0,10])
%     ylim([0,max(response)+(max(response)*0.1)]);
    
%     line([0,55],[threshold,threshold])
    
    % Go earlier
    earlier = 1;
    start_position = max_response_position(c);
    end_position = max_response_position(c);
    
    try
        while earlier
            if response(max_response_position(c)-earlier) > threshold
                start_position = max_response_position(c)-earlier;
                earlier = earlier+1;
            else
                earlier = 0;
            end
        end
    end
    later = 1;
    
    try
        while later
             if response(max_response_position(c)+later) > threshold
                end_position = max_response_position(c)+later;
                later = later+1;
            else
                later = 0;
             end
        end
    end
    
%     line([start_position,start_position],[0,10]);
%     line([end_position,end_position],[0,10]);
    
    lcs_field_width(c) = end_position-start_position;
    
end


% Place cells
pcs = [];
for s = 1:length(mData)
        
    if isempty(pcs)
        pcs = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on);
    else
        pcs = [pcs; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on)];
    end

    
end


pcs_field_width = [];

% Find field width for each cell at the first landmark position
for c = 1:size(pcs,1)
    response = pcs(c,:);
    [~,max_response_position(c)] = nanmax(response);
    threshold = nanmean(response);
    
    % Go earlier
    earlier = 1;
    start_position = max_response_position(c);
    end_position = max_response_position(c);
    
    try
        while earlier
            if response(max_response_position(c)-earlier) > threshold
                start_position = max_response_position(c)-earlier;
                earlier = earlier+1;
            else
                earlier = 0;
            end
        end
    end
    later = 1;
    
    try
        while later
             if response(max_response_position(c)+later) > threshold
                end_position = max_response_position(c)+later;
                later = later+1;
            else
                later = 0;
             end
        end
    end
    
%     line([start_position,start_position],[0,10]);
%     line([end_position,end_position],[0,10]);

    pcs_field_width(c) = end_position-start_position;
   
end

figure(1);
subplot(1,2,2);
landmark_cells_color = [239,35,60]/255;
place_cells_color = [51,51,51]/255;

hold on;
histogram(pcs_field_width,'BinEdges',0:2:50,'Normalization','probability','FaceColor',place_cells_color,'FaceAlpha',1)
hold on
histogram(lcs_field_width,'BinEdges',0:2:50,'Normalization','probability','FaceColor',landmark_cells_color,'FaceAlpha',0.8)

yticks([0,0.05,0.1,0.15,0.2])
yticklabels({0,5,10,15,20})
ylabel('% of cells');
ylim([0,0.18])
xticks([0,50])
xlim([0,50])
xlabel('Width (cm)')

title('Field width distribution');

set(gca,'FontSize',16);

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [1000,100,500,250])

[~ ,p_amp]= kstest2(amp_pcs,amp_lcs)
[~,p_width]= kstest2(lcs_field_width,pcs_field_width)



function [mData,data] = rsc_jitter()
dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/Nora/jittering/';

% Recordings to include
data.sessionIDs = { 'm1551-20230308-03', ... # Landmarks removed after 41 laps. Ran 96 laps total
                    'm1551-20230309-01', ... # Landmarks removed after 41 laps. Ran 136 laps total
                    'm1553-20230308-01', ... # Landmarks removed after 30 laps. Ran 44 laps total
                    'm1553-20230309-01', ...
                    'm1554-20230308-01',...
                    'm1554-20230309-01'
                    };
                
                
% % Load all sessions

for f = 1:length(data.sessionIDs)
    
    load(fullfile(dirct,data.sessionIDs{f},'sData.mat'));
    mData(f).sData = sData;
end



% Initiate some variables
% data.n_cells_pos1 = [];
data.n_lcs_pos1 = 0;
data.n_lcs_pos2 = 0;
data.n_lcs_pos3 = 0;
data.n_lcs_pos4 = 0;

max_z_value = 3;

% Create rastermaps for each session
for s = 1:length(mData)
    
    mData(s).sData.imdata.roiSignals(2).('dffSubtractedNpil') = mData(s).sData.imdata.roiSignals(2).dff;
    
    mData(s).sData.rmaps = split_rmaps_jitter(mData(s).sData);

    % identify landmark and place cells in on trials 
 
    params.first_landmark_bins = 15:29; params.second_landmark_bins = 48:62;
    mData(s).rmaps_first = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_firstpair,params);  
    params.first_landmark_bins = 25:39; params.second_landmark_bins = 58:72;    
    mData(s).rmaps_second = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_secondpair,params); 
    params.first_landmark_bins = 35:48; params.second_landmark_bins = 68:82;   
    mData(s).rmaps_third = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_thirdpair,params); 
    params.first_landmark_bins = 45:59; params.second_landmark_bins = 78:92;  
    mData(s).rmaps_fourth = sb.classify.landmarkCellsByShuffling(mData(s).sData.rmaps.deconv_fourthpair,params);  

   data.n_lcs_pos1 = data.n_lcs_pos1 +length(mData(s).rmaps_first);
   data.n_lcs_pos2 = data.n_lcs_pos2 + length(mData(s).rmaps_second);
   data.n_lcs_pos3 = data.n_lcs_pos3 + length(mData(s).rmaps_third);
   data.n_lcs_pos4 = data.n_lcs_pos4 + length(mData(s).rmaps_fourth);
    
   lcs_in_12 = intersect(mData(s).rmaps_first,mData(s).rmaps_second);
   lcs_in_123 =  intersect(lcs_in_12,mData(s).rmaps_third);
   data.lcs_in_all{s} =  intersect(lcs_in_123,mData(s).rmaps_fourth);
   data.n_lcs_in_all(s) = length(data.lcs_in_all{s});
   data.prct_lcs_in_all(s) = length(data.lcs_in_all{s})/size(mData(s).sData.rmaps.deconv_firstpair,3);

end



figure()

scatter(ones(1,6)',100*data.prct_lcs_in_all',50,'k')
hold on
scatter(1,100*nanmean(data.prct_lcs_in_all'),70,'k','filled')
errorbar(1,100*nanmean(data.prct_lcs_in_all'),nanstd(100*data.prct_lcs_in_all')/sqrt(length(data.prct_lcs_in_all)),'k','LineWidth',1.5)
ylabel('Object resposinve cells in all 4 (%)')
set(gca,'FontName','Arial', 'FontSize',14)
xticks([])
yticks([0:10:50])


amap_object_responsive_in_all_pos1 = [];
amap_object_responsive_in_all_pos2 = [];
amap_object_responsive_in_all_pos3 = [];
amap_object_responsive_in_all_pos4 = [];
for s = 1:length(mData)

     amap_object_responsive_in_all_pos1 = cat(1,amap_object_responsive_in_all_pos1,sb.generate.averagedTuningFromRasterMaps(mData(s).sData.rmaps.deconv_firstpair(:,:,data.lcs_in_all{s})));
     amap_object_responsive_in_all_pos2 = cat(1,amap_object_responsive_in_all_pos2,sb.generate.averagedTuningFromRasterMaps(mData(s).sData.rmaps.deconv_secondpair(:,:,data.lcs_in_all{s})));
     amap_object_responsive_in_all_pos3 = cat(1,amap_object_responsive_in_all_pos3,sb.generate.averagedTuningFromRasterMaps(mData(s).sData.rmaps.deconv_thirdpair(:,:,data.lcs_in_all{s})));
     amap_object_responsive_in_all_pos4 = cat(1,amap_object_responsive_in_all_pos4,sb.generate.averagedTuningFromRasterMaps(mData(s).sData.rmaps.deconv_fourthpair(:,:,data.lcs_in_all{s})));
    
end

[amap_object_responsive_in_all_pos1,sortIdx_lcs]  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos1 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos1  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos1 ,'normalize_type','zscore');

[amap_object_responsive_in_all_pos2,sortIdx_lcs]  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos2 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos2  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos2 ,'normalize_type','zscore');

[amap_object_responsive_in_all_pos3,sortIdx_lcs]  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos3 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos3  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos3 ,'normalize_type','zscore');

[amap_object_responsive_in_all_pos4,sortIdx_lcs]  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos4 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos4  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos4 ,'normalize_type','zscore');

figure()
subplot(1,4,1)
imagesc(amap_object_responsive_in_all_pos1,[0,max_z_value])
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([22,22],[0,2000],'Color','w');
line([55,55],[0,2000],'Color','w');
ylabel('Cell #')


subplot(1,4,2)
imagesc(amap_object_responsive_in_all_pos2,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([32,32],[0,2000],'Color','w');
line([65,65],[0,2000],'Color','w');

subplot(1,4,3)
imagesc(amap_object_responsive_in_all_pos3,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([42,42],[0,2000],'Color','w');
line([75,75],[0,2000],'Color','w');

subplot(1,4,4)
imagesc(amap_object_responsive_in_all_pos4,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([52,52],[0,2000],'Color','w');
line([85,85],[0,2000],'Color','w');

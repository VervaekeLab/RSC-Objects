
% PAPER 1 FIG 2
% Plot a touch cell / landmark cell that doesn't last, and a lasting
% touch-cell / trace cell

figure(2);
dot_size = 3;
line_width = 1;
max_z_value = 4; 

% Motor on marks
landmark_onsets = mData(5).sData.landmarks.motor.on(:,1);

% Non-lasting touch cell
subplot(1,2,1)
rmap = [mData(5).rmaps.lcs.deconv_motor_on(:,:,8);mData(5).rmaps.lcs.deconv_motor_off(:,:,8)];
rmap = normalize(rmap);
% rmap = zScoreNormalize(rmap,'all');

imagesc(rmap,[0,max_z_value])
colormap(flipud(colormap('gray')));
hold on;
line([32,32],[0,size(mData(5).rmaps.lcs.deconv_motor_on(:,:,8),1)],'LineWidth',1,'Color','r')
line([72,72],[0,1000],'LineWidth',1,'Color','r')
line([0,105],[nansum(landmark_onsets),nansum(landmark_onsets)],'LineWidth',line_width,'Color','k')

xticks([]);
yticks([1,size(mData(5).rmaps.lcs.deconv_motor_on(:,:,8),1),size(rmap,1)])
ylabel('Lap #');
set(gca,'FontSize',16)

% Lasting touch cell
subplot(1,2,2)
rmap = [mData(5).rmaps.tcs.deconv_motor_on(:,:,4);mData(5).rmaps.tcs.deconv_motor_off(:,:,4)];
% rmap = zScoreNormalize(rmap,'all');
rmap = normalize(rmap);
imagesc(rmap,[0,max_z_value])
colormap(flipud(colormap('gray')));
hold on;

line([32,32],[0,size(mData(5).rmaps.lcs.deconv_motor_on(:,:,4),1)],'LineWidth',1,'Color','r')
line([72,72],[0,1000],'LineWidth',1,'Color','r')
line([0,105],[nansum(landmark_onsets),nansum(landmark_onsets)],'LineWidth',line_width,'Color','k')

xticks([]);
yticks([])

set(gcf,'renderer', 'painters', 'Position', [2000 200 500 300])
function alignmentAndBehavior(sData)
%% Plot the alignment vector together with running speed and postion of animal
% Useful to look for movement problems in image data.

figure(1);clf;


max_length = min([length(sData.imdata.alignment.nonrigid_corrections(:,1)),length(sData.behavior.wheelPosDs)]);

subplot(4,1,1);
plot(sData.imdata.alignment.nonrigid_corrections(:,1),'Color','r')
hold on
plot(sData.imdata.alignment.nonrigid_corrections(:,2),'Color','r')
title('Alignment');
xlim([1,max_length]);

subplot(4,1,2);
plot(sData.behavior.wheelPosDs,'Color','k','LineWidth',0.5);

max_movements_original = max(abs(sData.imdata.alignment.nonrigid_corrections(:,:))');
max_movements = max_movements_original;
max_movements(max_movements<10) = 0;
hold on;
plot(normalize(max_movements,'range',[0,157]))

xlim([1,max_length]);
title('Position')

subplot(4,1,3);
plot(sData.behavior.runSpeedDs,'Color','k','LineWidth',2);
title('Running speed')
xlim([1,max_length]);

subplot(4,1,4);

rmap = sb.generate.rasterMaps(sData.behavior.wheelPosDsBinned,sData.behavior.wheelLapDsBinned,max_movements_original);
imagesc(normalize(rmap,'zscore'),[0,5]);
title('Rastermap of max shifts vs position');

end
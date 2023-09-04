function roiImageFromMultipleSessions(mData,roi)
% Plot the roi image for a given roi number across all session in mData
% given as input.


n_plots = length(mData);
figure(23121);
clf;

for p = 1:n_plots
    subplot(1,n_plots,p);
    
    imagesc(mData(p).sData.imdata.roiArray(roi).enhancedImage);
    colormap('gray');
    title(sprintf('%s',strrep(mData(p).sData.sessionInfo.sessionID(1:17),'_','-')),'FontSize',16)
end





end
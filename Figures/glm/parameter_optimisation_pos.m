N_P = 1:8;%
beta = [0 0.01 0.1 1 10];
[f,n]     = ndgrid(N_P,beta);
Z         = [f(:) n(:)];

                          
for i = 1:length(data.sessionIDs)
%     load(fullfile(dirct,data.sessionIDs{i},'sData.mat'))
    for j = 1:length(Z)
        if exist(fullfile(dirct_optim,data.sessionIDs{i},strcat('pos_hyper_par_',num2str(j),'.mat')),'file')
            load(fullfile(dirct_optim,data.sessionIDs{i},strcat('pos_hyper_par_',num2str(j),'.mat')))
            for m = 1:length(axon)
               LLH{i}(m,j,:) = axon(m).testFit{1,j}(:,3);
            end
            
        else
            LLH{i}(m,j,:) = NaN(5,1);
        end
        
    end

    for l = 1:length(N_P)
        for j =1:length(beta)
            idx = find(Z(:,1)==N_P(l));
            idxS = find(Z(idx,2) == beta(j));
            LLH_mean(i,l,j) = nanmean(nanmean(LLH{i}(:,idx(idxS),:),3),1);

        end
    end
    
    [X,Y]= meshgrid(N_P,beta);
    fig = figure();
    surf(X',Y',reshape(LLH_mean(i,:,:),8,5), 'edgecolor', 'none')
    xlabel('Position bin (cm)')
    xticks([1:8])
    xticklabels([1.5:1.5:1.5*8])
    ylabel('smoothing')
    zlabel('LLH-pos')
    colormap('jet')
    set(gca,'yscale','log')
    
    
end

fig = figure();
surf(X',Y',reshape(nanmean(LLH_mean,1),8,5), 'edgecolor', 'none','FaceAlpha',0.4)
xlabel('Position bin (cm)')
xticks([1:8])
xticklabels([1.5:1.5:1.5*8])
ylabel('smoothing')
zlabel('LLH-pos')
colormap('jet')
set(gca,'yscale','log')


ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];


LLH_mean_all = reshape(nanmean(LLH_mean,1),8,5);
[maxSmooth,binIdx] = max(LLH_mean_all);
[maxBinSize, smoothIdx]  = max(maxSmooth);
maxbinIdx = binIdx(smoothIdx);
hold on
scatter3(N_P(maxbinIdx),beta(smoothIdx),LLH_mean_all(maxbinIdx,smoothIdx),50,'k','filled')
    


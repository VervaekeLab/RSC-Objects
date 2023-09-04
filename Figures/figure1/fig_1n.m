
% this analysis was run on a super computer to speed up analysis time
% parameters were optimised, the optisation analysis can be viewed here:
% /Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/SVM/paramOptim
% optimised parameters were : time binning  = 0.1935 (6 time samples)
%                             learning rate = 0.1
% 
% epoches = 300;
Epoches = 300;

runSVM = false;

if runSVM
    % SVM_pos1pos2.m % use this function, needs still to include loop over sData! 
else
    dirct_SVM = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/SVM/pos1pos2';
    for s = 1:length(data.sessionIDs)
        load(fullfile(dirct_SVM,data.sessionIDs{s},'d_data.mat'));
        d_data(s) = d_data_temp(1);
    end
end
    



crossvalidationFolds = 5;
mean_total_accuracy = [];
mean_class   = []; 
accuracy_class_fold = [];
total_accuracy_fold = [];

for s = 1:length(data.sessionIDs)
    for m = 1:5
        [accuracy_class_fold(s,m,:),total_accuracy_fold(s,m)]= classificationAccuracy(d_data(s).iter{1, m}.per_outlabel_test{Epoches}  , d_data(s).iter{1, m}.realPos_test,[-1 1]);
    end
    mean_total_accuracy(s)   = nanmean(total_accuracy_fold(s,:),2);
    mean_class(s,:)          = nanmean(accuracy_class_fold(s,:,:),2);
end


fig = figure();

scatter(ones(1,length(data.sessionIDs)),mean_total_accuracy,50,'k');hold on
hold on
scatter(1,nanmean(mean_total_accuracy,2),50,'MarkerFaceColor','k','MarkerEdgeColor','none'); hold on
errorbar(1,nanmean(mean_total_accuracy,2),nanstd(mean_total_accuracy',1)'/sqrt(length(data.sessionIDs)),'LineStyle','none','LineWidth',1,'Color','k')

ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
xlim([0.5 1.5])
ylim([0.4 0.8])
yticks([0.4:0.1:0.8])
yline(0.5)
yticklabels([40:10:80])
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Prediction accuracy (%)')
xticks([1])
xticklabels({'First landmark vs second landmark'})
xtickangle(30)

    


function [classAcc,total_accuracy]= classificationAccuracy(predClass, actualClass,classes)

classAcc = zeros(length(classes),1);

for i = 1: length(classes)
    idx                   = actualClass == classes(i);
    predClasstemp         = predClass(idx);
    corrPred(i)           = numel(find(predClasstemp == classes(i)));
    classAcc(i)           = corrPred(i)/length(find(idx==1));
    
end

total_accuracy= sum(corrPred)/length(actualClass);

end



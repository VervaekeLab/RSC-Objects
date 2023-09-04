
% this analysis was run on a super computer to speed up analysis time
% parameters were optimised, the optisation analysis can be viewed here:
% /Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/SVM/paramOptim
% optimised parameters were : time binning  = 0.1935 (6 time samples)
%                             learning rate = 0.1
% 


runSVM = false;

if runSVM
    % run sb.svm_pos1pos2_type1 for all sData
    % run sb.svm_pos1pos2_type1 for all sData
    % run sb.svm_type1type2 for all sData
else
    dirct_SVM = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/SVM/two_objects/';
    predType = {'pos1pos2_type1_only_lcsA','pos1pos2_type2_only_lcsB','type1type2_only_lcsAB_equalNoAB'};

    for s = 1:length(data.sessionIDs)
        load(fullfile(dirct_SVM,'pos1pos2_type1_only_lcsA',data.sessionIDs{s},'d_data.mat'));
        d_data_pos1pos2_type1(s) = d_data_temp(1);
    end

    for s = 1:length(data.sessionIDs)
        load(fullfile(dirct_SVM,'pos1pos2_type2_only_lcsB',data.sessionIDs{s},'d_data.mat'));
        d_data_pos1pos2_type2(s) = d_data_temp(1);
    end
    
    for s = 1:length(data.sessionIDs)
        load(fullfile(dirct_SVM,'type1type2_all_lcs',data.sessionIDs{s},'d_data.mat'));
        d_data_type1type2(s) = d_data_temp(1);
    end
    
end
   

crossvalidationFolds = 5;
mean_total_accuracy = [];
mean_class   = []; 
accuracy_class_fold = [];
total_accuracy_fold = [];
xticklabels({'First landmark vs second landmark'})
xtickangle(30)

    
for s = 1:length(data.sessionIDs)
    for m = 1:5
        [~,total_accuracy(s,m)]= classificationAccuracy(d_data_pos1pos2_type1(s).iter{1, m}.per_outlabel_test{500},d_data_pos1pos2_type1(s).iter{m}.realPos_test,[-1 1]);
    end
    accuracy(2,s)   = nanmean(total_accuracy(s,:),2);
end


for s = 1:length(data.sessionIDs)
    for m = 1:5
        [~,total_accuracy(s,m)]= classificationAccuracy(d_data_pos1pos2_type2(s).iter{1, m}.per_outlabel_test{500},d_data_pos1pos2_type2(s).iter{m}.realPos_test,[-1 1]);
    end
    accuracy(1,s)   = nanmean(total_accuracy(s,:),2);
end

for s = 1:length(data.sessionIDs)
    for m = 1:5
        [~,total_accuracy(s,m)]= classificationAccuracy(d_data_type1type2(s).iter{1, m}.per_outlabel_test{500},d_data_type1type2(s).iter{m}.realPos_test,[-1 1]);
    end
    accuracy(3,s)   = nanmean(total_accuracy(s,:),2);
end


% for s = 1:length(data.sessionIDs)
%     no_sandpaper(s) = length(mData(s).lcs_sandpaper);
%     no_feltpad(s) = length(mData(s).lcs_feltpad);
% end
% 
% accuracy(2,no_sandpaper<5) = NaN;
% accuracy(1,no_feltpad<5) = NaN;
% accuracy(3,(no_feltpad+no_sandpaper)<5) = NaN;

fig = figure();
for type = 1:length(predType)
    scatter(type*ones(1,length(data.sessionIDs)),accuracy(type,:),70,'k');hold on
    scatter(type,nanmean(accuracy(type,:)),70,'MarkerFaceColor','k','MarkerEdgeColor','none'); hold on
    errorbar(type,nanmean(accuracy(type,:)),nanstd(accuracy(type,:))/sqrt(sum(~isnan(accuracy(type,:)))),'LineStyle','none','LineWidth',1,'Color','k')
end


ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

xticks([1 2 3])
xlim([0.5 3.5])
ylim([0.4 0.8])
yline(0.5)

xticklabels({'Felt pad positions','Sandpaper positions','Feltpad vs sandpaper'});
xtickangle(45)
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Prediction accuracy')
% 
pval(1) = signrank(0.5*ones(9,1),accuracy(1,:)');
pval(2) = signrank(0.5*ones(9,1),accuracy(2,:)');
pval(3) = signrank(0.5*ones(9,1),accuracy(3,:)');

yline(0.5)






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

function varargout=sigstar(groups,stats,nosort)
    % sigstar - Add significance stars to bar charts, boxplots, line charts, etc,
    %
    % H = sigstar(groups,stats,nsort)
    %
    % Purpose
    % Add stars and lines highlighting significant differences between pairs of groups. 
    % The user specifies the groups and associated p-values. The function handles much of 
    % the placement and drawing of the highlighting. Stars are drawn according to:
    %   * represents p<=0.05
    %  ** represents p<=1E-2
    % *** represents p<=1E-3
    %
    %
    % Inputs
    % groups - a cell array defining the pairs of groups to compare. Groups defined 
    %          either as pairs of scalars indicating locations along the X axis or as 
    %          strings corresponding to X-tick labels. Groups can be a mixture of both 
    %          definition types.
    % stats -  a vector of p-values the same length as groups. If empty or missing it's 
    %          assumed to be a vector of 0.05s the same length as groups. Nans are treated
    %          as indicating non-significance.
    % nsort -  optional, 0 by default. If 1, then significance markers are plotted in 
    %          the order found in groups. If 0, then they're sorted by the length of the 
    %          bar.
    %
    % Outputs
    % H - optionally return handles for significance highlights. Each row is a different
    %     highlight bar. The first column is the line. The second column is the text (stars).
    %     
    %
    % Examples
    % 1. 
    % bar([5,2,1.5])
    % sigstar({[1,2], [1,3]})
    %
    % 2. 
    % bar([5,2,1.5])
    % sigstar({[2,3],[1,2], [1,3]},[nan,0.05,0.05])
    %
    % 3.  **DOESN'T WORK IN 2014b**
    % R=randn(30,2);
    % R(:,1)=R(:,1)+3;
    % boxplot(R)
    % set(gca,'XTick',1:2,'XTickLabel',{'A','B'})
    % H=sigstar({{'A','B'}},0.01);
    % ylim([-3,6.5])
    % set(H,'color','r')
    %
    % 4. Note the difference in the order with which we define the groups in the 
    %    following two cases. 
    % x=[1,2,3,2,1];
    % subplot(1,2,1)
    % bar(x)
    % sigstar({[1,2], [2,3], [4,5]})
    % subplot(1,2,2)
    % bar(x)
    % sigstar({[2,3],[1,2], [4,5]})
    %
    % ALSO SEE: demo_sigstar
    %
    % KNOWN ISSUES:
    % 1. Algorithm for identifying whether significance bar will overlap with 
    %    existing plot elements may not work in some cases (see line 277)
    % 2. Bars may not look good on exported graphics with small page sizes.
    %    Simply increasing the width and height of the graph with the 
    %    PaperPosition property of the current figure should fix things.
    %
    % Rob Campbell - CSHL 2013
    %Input argument error checking
    %If the user entered just one group pair and forgot to wrap it in a cell array 
    %then we'll go easy on them and wrap it here rather then generate an error
    if ~iscell(groups) & length(groups)==2
        groups={groups};
    end
    if nargin<2 
        stats=repmat(0.05,1,length(groups));
    end
    if isempty(stats)
        stats=repmat(0.05,1,length(groups));
    end
    if nargin<3
        nosort=0;
    end
    %Check the inputs are of the right sort
    if ~iscell(groups)
        error('groups must be a cell array')
    end
    if ~isvector(stats)
        error('stats must be a vector')
    end
    if length(stats)~=length(groups)
        error('groups and stats must be the same length')
    end
    %Each member of the cell array groups may be one of three things:
    %1. A pair of indices.
    %2. A pair of strings (in cell array) referring to X-Tick labels
    %3. A cell array containing one index and one string
    %
    % For our function to run, we will need to convert all of these into pairs of
    % indices. Here we loop through groups and do this. 
    xlocs=nan(length(groups),2); %matrix that will store the indices 
    xtl=get(gca,'XTickLabel');  
    for ii=1:length(groups)
        grp=groups{ii};
        if isnumeric(grp)
            xlocs(ii,:)=grp; %Just store the indices if they're the right format already
        elseif iscell(grp) %Handle string pairs or string/index pairs
            if isstr(grp{1})
                a=strmatch(grp{1},xtl);
            elseif isnumeric(grp{1})
                a=grp{1};
            end
            if isstr(grp{2})
                b=strmatch(grp{2},xtl);
            elseif isnumeric(grp{2})
                b=grp{2};
            end
            xlocs(ii,:)=[a,b];
        end
        %Ensure that the first column is always smaller number than the second
        xlocs(ii,:)=sort(xlocs(ii,:));
    end
    %If there are any NaNs we have messed up. 
    if any(isnan(xlocs(:)))
        error('Some groups were not found')
    end
    %Optionally sort sig bars from shortest to longest so we plot the shorter ones first
    %in the loop below. Usually this will result in the neatest plot. If we waned to 
    %optimise the order the sig bars are plotted to produce the neatest plot, then this 
    %is where we'd do it. Not really worth the effort, though, as few plots are complicated
    %enough to need this and the user can define the order very easily at the command line. 
    if ~nosort
        [~,ind]=sort(xlocs(:,2)-xlocs(:,1),'ascend');
        xlocs=xlocs(ind,:);groups=groups(ind);
        stats=stats(ind);
    end
    %-----------------------------------------------------
    %Add the sig bar lines and asterisks 
    holdstate=ishold;
    hold on
    H=ones(length(groups),2); %The handles will be stored here
    y=ylim;
    yd=myRange(y)*0.05; %separate sig bars vertically by 5% 
    for ii=1:length(groups)
        thisY=findMinY(xlocs(ii,:))+yd;
        H(ii,:)=makeSignificanceBar(xlocs(ii,:),thisY,stats(ii));
    end
    %-----------------------------------------------------
    %Now we can add the little downward ticks on the ends of each line. We are
    %being extra cautious and leaving this it to the end just in case the y limits
    %of the graph have changed as we add the highlights. The ticks are set as a
    %proportion of the y axis range and we want them all to be the same the same
    %for all bars.
    yd=myRange(ylim)*0.03; %Ticks are 1% of the y axis range
    for ii=1:length(groups)
        y=get(H(ii,1),'YData');
        y(1)=y(1)-yd;
        y(4)=y(4)-yd;   
        set(H(ii,1),'YData',y)
    end
    %Be neat and return hold state to whatever it was before we started
    if ~holdstate
        hold off
    elseif holdstate
        hold on
    end
    %Optionally return the handles to the plotted significance bars (first column of H)
    %and asterisks (second column of H).
    if nargout>0
        varargout{1}=H;
    end
end %close sigstar

function H=makeSignificanceBar(x,y,p)
    %makeSignificanceBar produces the bar and defines how many asterisks we get for a 
    %given p-value
    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='n.s.';
    end
            
    x=repmat(x,2,1);
    y=repmat(y,4,1);
    H(1)=plot(x(:),y,'-k','LineWidth',1.5,'Tag','sigstar_bar');
    %Increase offset between line and text if we will print "n.s."
    %instead of a star. 
%     if ~isnan(p)
%         offset=0.005;
%     else
        offset=0.02;
%     end
    starY=mean(y)+myRange(ylim)*offset;
    H(2)=text(mean(x(:)),double(starY),stars,...
        'HorizontalAlignment','Center',...
        'BackGroundColor','none',...
        'Tag','sigstar_stars','FontSize',16);
    Y=ylim;
    if Y(2)<starY
        ylim([Y(1),starY+myRange(Y)*0.05])
    end
end %close makeSignificanceBar

function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');
    axis(gca,'tight')
    
    %increase range of x values by 0.1 to ensure correct y max is used
    x(1)=x(1)-0.1;
    x(2)=x(2)+0.1;
    
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis
    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);
    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)
end %close findMinY

function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end %close myRange


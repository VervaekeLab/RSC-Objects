


clear all
close all
% this function needs :
% sData struct
% lcs: a vector of object location cell indices
% tcs: a vector of trace cell indices

% here we use a time bin size of 6 and a learning rate of 0.1, these were
% optimised

timeBinSize = 6;
learningRate = 0.1;

file  = '';
save_dirct = '';

if ~exist(save_dirct, 'dir'); mkdir(save_dirct); end

lcs= unique([lcs,tcs]);

d_data = [];

deconvSignal=[];
for i = 1:size(sData.imdata.roiSignals(2).dff,1)
    deconvSignal(i,:) = double(normalize(smoothdata(sData.imdata.roiSignals(2).deconv(i,:),5),'range',[0 1]));
end

position        = sData.behavior.wheelPosDs;
runSpeed        = sData.behavior.runSpeedDs;
lap     = sData.behavior.wheelLapDsBinned;
% Shift lap to start on 1, not -1 or 0.
lap = lap + (abs(min(lap))*2);

% Find laps where the landmark motor is present (on) or omitted (off)
laps_landmark_on  = find(sData.landmarks.motor.on(:,1) == 1);
on_laps  = find(ismember(lap,laps_landmark_on));

position        = position(on_laps );
deconvSignal    = deconvSignal(lcs',on_laps);
runSpeed        = runSpeed(on_laps);
cellIdx         = lcs;


deconvSignalBinned  = NaN(size(deconvSignal,1),length(timeBins)-1);
runSpeedBinned      = NaN(length(timeBins)-1,1);
positionBinned      = NaN(length(timeBins)-1,1);
timeBins = 1:timeBinSize:length(position);
for t = 1:length(timeBins)-1
    positionBinned(t)        = nanmean(position(timeBins(t):timeBins(t+1)-1));
    deconvSignalBinned(:,t)  = nanmean(deconvSignal(:,timeBins(t):timeBins(t+1)-1),2);
    runSpeedBinned(t)        = nanmean(runSpeed(timeBins(t):timeBins(t+1)-1));
end

% delete all time points in whhhichh the aanimal runs slower than 1 cm/s
positionBinned(runSpeedBinned<1)        = [];
deconvSignalBinned(:,runSpeedBinned<1)  = [];
runSpeedBinned(runSpeedBinned<1)        = [];

noCells     = size(deconvSignalBinned,1);
posBinned   = ones(length(positionBinned),1);

idxFirst    = find(positionBinned >= 34 & positionBinned <= 65);
idxSecond   = find(positionBinned >= 94 & positionBinned <= 125);

minLength   = min(length(idxFirst),length(idxSecond));
idxFirst    = randsample(idxFirst,minLength);
idxSecond   = randsample(idxSecond,minLength);

posBinnedTemp(idxFirst)     =  1;
posBinnedTemp(idxSecond)    = -1;
bins = 1:2;

% remove all indices that are not around landmarks
idxNotUsed = setdiff(1:length(posBinnedTemp),[idxFirst(:);idxSecond(:)]);
deconvSignalBinned(:,idxNotUsed) = [];
posBinnedTemp(idxNotUsed) = [];

% divide the data up into 5*num_folds pieces
numFolds = 5;
sections = numFolds*5;
edges = round(linspace(1,size(deconvSignalBinned,2)+1,sections+1));

rmse_train = zeros(1,numFolds);
rmse_test = zeros(1,numFolds);

for k = 1:numFolds
    testIdx  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1];
    
    trainIdx = setdiff(1:size(deconvSignalBinned,2),testIdx);
    
    d_data.iter{k}.deconvTrain = deconvSignalBinned(:,trainIdx);
    d_data.iter{k}.deconvTest  = deconvSignalBinned(:, testIdx);
    %
    %
    %% find cells without signal and delete
    d_data.iter{k}.cellIdx = cellIdx;
    d_data.iter{k}.realPos_train = posBinnedTemp(trainIdx); YTrain =  posBinnedTemp(trainIdx)';
    d_data.iter{k}.realPos_test  = posBinnedTemp(testIdx);  YTest =  posBinnedTemp(testIdx)';
    XTrain = d_data.iter{k}.deconvTrain;
    XTest  = d_data.iter{k}.deconvTest;
    
    
    [d_data.iter{k}.per,d_data.iter{k}.per_outlabel_test, d_data.iter{k}.per_outlabel_train] = cuesSVM(XTrain, YTrain, XTest, YTest,learningRate(m));
    
end
d_data_temp = d_data;
save(fullfile(dirct, strcat('/d_data','.mat')), 'd_data_temp','-v7.3');





function [per,per_outlabel_test,per_outlabel_train] = cuesSVM(XTrain, YTrain, XTest, YTest,sigvals)  
%% Perceptron Algorithm
%% Initialization
count=0; temp=0;

train   = XTrain';
test    = XTest';
labeltr = YTrain';
labelte = YTest';

summation=0;
N1_max=length(YTrain);
alpha=zeros(N1_max,1);
iter=1; y_hat=0;y_pred=0;
% sigvals = 0.1;
for N1=N1_max%2000:1000:N1_max
    for l=1:length(sigvals)
        sig = sigvals(l);
        alpha=zeros(N1,1);
        N2=length(YTest);%abs(0.04*N1); %Number of data points from testing set
        %% Learning
        for n=1:300 %Number of passes over the data
            tic
            for j=1:N1
                for i=1:N1
                    y_hat=y_hat+alpha(i,1)*labeltr(i,1)*gaussian_kernel(train(i,:),train(j,:),sig);
                end
                if(sign(y_hat)~=labeltr(j))
                    
                    alpha(j)=alpha(j)+1;
                end
                y_hat=0;
            end
            t=toc;
            
            %% Prediction test data
            for j=1:N2
                for i=1:N1
                    y_pred=y_pred+alpha(i,1)*labeltr(i)*gaussian_kernel(train(i,:),test(j,:),sig);
                end
                
                if(sign(y_pred)>=0)
                    outlabelTest(j,1)=1;
                else
                    outlabelTest(j,1)=-1;
                end
                y_pred=0;
            end
            
            %% Prediction test data
            for j=1:N1
                for i=1:N1
                    y_pred=y_pred+alpha(i,1)*labeltr(i)*gaussian_kernel(train(i,:),train(j,:),sig);
                end
                
                if(sign(y_pred)>=0)
                    outlabelTrain(j,1)=1;
                else
                    outlabelTrain(j,1)=-1;
                end
                y_pred=0;
            end
            
            
            %% Accuracy
            A=(outlabelTest==labelte(1:N2));
            count=sum(A(:));
            per{l,iter}(1)=N1;
            if n==1
                per{l,iter}(2)=t;
            else
                per{l,iter}(2)=per{l,iter-1}(2)+t;
            end
            per{l,iter}(3)=sigvals(l);
            per{l,iter}(4)=n;
            per{l,iter}(5)=count/N2*100;
            per_outlabel_test{l,iter} = outlabelTest;
            per_outlabel_train{l,iter} = outlabelTrain;
            per{l,iter}(6) = count/length(A);
            iter=iter+1;
            count=0;
        end
    end
    
end

end
    
function [dotp]=gaussian_kernel(x,y,sig)
if nargin<3
sig=0.2;
end

dotp=exp(-(norm(x-y)^2)/(2*sig^2));
end
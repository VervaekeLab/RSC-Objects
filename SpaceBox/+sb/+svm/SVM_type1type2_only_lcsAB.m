


clear all
close all

% this function needs :
% sData struct
% lcsA: a vector of object A cell indices
% lcsB: a vector of object B cell indices

% here we use a time bin size of 6 and a learning rate of 0.1, these were
% optimised

timeBinSize = 6;
learningRate = 0.1;

file  = '';
save_dirct = '';

if ~exist(save_dirct, 'dir'); mkdir(save_dirct); end

lcs= unique([lcsA,lcsB]);

d_data = [];
deconvSignal=[];
for i = 1:size(sData.imdata.roiSignals(2).dff,1)
    deconvSignal(i,:) = double(normalize(smoothdata(sData.imdata.roiSignals(2).deconv(i,:),5),'range',[0 1]));
end

position        = sData.behavior.wheelPosDs;
runSpeed        = sData.behavior.runSpeedDs;

% make sure all variables are of same length
minLength       = min([length(position), size(deconvSignal,2), length(runSpeed)]);
position        = position(1:minLength);
deconvSignal    = deconvSignal(lcs',1:minLength);
runSpeed        = runSpeed(1:minLength);
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
posBinned            = ones(length(positionBinned),1);

idxFirst = [find(positionBinned >= 42-8& positionBinned <= 47+8),find(positionBinned >= 122-8 & positionBinned <= 128+8)];
idxSecond = [find(positionBinned >= 63-8& positionBinned <= 68+8),find(positionBinned >= 100-8 & positionBinned <= 106+8)];

minLength = min(length(idxFirst),length(idxSecond));
idxFirst = randsample(idxFirst,minLength);
idxSecond = randsample(idxSecond,minLength);


minLength = min(length(idxFirst),length(idxSecond));
idxFirst = randsample(idxFirst,minLength);
idxSecond = randsample(idxSecond,minLength);

posBinned(idxFirst) =  1;
posBinned(idxSecond) = -1;
bins = 1:2;

% remove all indices that are not around landmarks
idxNotUsed = setdiff(1:length(posBinned),[idxFirst(:);idxSecond(:)]);
deconvSignalBinned(:,idxNotUsed) = [];
posBinned(idxNotUsed) = [];

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
    
    d_d_data.iter{k}.deconvTrain = deconvSignalBinned(:,trainIdx);
    d_d_data.iter{k}.deconvTest  = deconvSignalBinned(:, testIdx);

    
    %% find cells without signal and delete
    d_d_data.iter{k}.cellIdx = cellIdx;
    d_d_data.iter{k}.realPos_train = posBinned(trainIdx); YTrain =  posBinned(trainIdx)';
    d_d_data.iter{k}.realPos_test  = posBinned(testIdx);  YTest =  posBinned(testIdx)';
    XTrain = d_d_data.iter{k}.deconvTrain;
    XTest  = d_d_data.iter{k}.deconvTest;
    
    
    [d_d_data.iter{k}.per,d_d_data.iter{k}.per_outlabel_test, d_d_data.iter{k}.per_outlabel_train] = cuesSVM(XTrain, YTrain, XTest, YTest,learningRate(m));
    
end
d_data_temp = d_d_data;

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
N1_max=length(YTrain);%2000; % Number of data points from training set
alpha=zeros(N1_max,1);
iter=1; y_hat=0;y_pred=0;
% sigvals = 0.1;%[0.01 0.1 0.2 0.3];
for N1=N1_max%2000:1000:N1_max
    for l=1:length(sigvals)
        sig = sigvals(l);
        alpha=zeros(N1,1);
        N2=length(YTest);%abs(0.04*N1); %Number of data points from testing set
        %% Learning
        for n=1:500 %Number of passes over the data
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
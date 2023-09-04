%% Clear the workspace and load the data

clear all; close all; clc
%% 
% to find N_T (time binning) and beta, we ran this script with all
% combinations of:
% N_T = [0.05,0.1:0.2:1,1]
% beta = [0 0.01 0.1 1 10]
% and found the values that lead to maximal LLH values with paramter_optimisation_time script

% results are saved in neuron


file = '';
load_dirct = '';
save_dirct = '';

dirct =  fullfile(save_dirct,file);
if ~exist(dirct, 'dir'); mkdir(dirct); end


load(fullfile(load_dirct,file,'/sData.mat'));

% binning parameters: 
N_P = 1;
beta = 1;

  
deconv  = double(sData.imdata.roiSignals(2).deconv);
pos     = sData.behavior.wheelPosDsBinned;
lap     = sData.behavior.wheelLapDsBinned;
speed   = sData.behavior.runSpeedDs;

deconv  = deconv(find(sum(deconv,2)),:);

numNeuron = length(lcs);
% Remove all bins where running speed is < 1 cm.)
bins_to_include = find([speed<1] == 0);
deconv          = deconv(lcs,bins_to_include);
pos             = pos(bins_to_include);
lap             = lap(bins_to_include);

% Add non-zero noise to deconv
deconv = deconv + randi(10,size(deconv,1),size(deconv,2))*0.0000000000000001;

% Smooth the deconv signal (Mao does this in his code, but don't write so in the paper)
deconv = smoothdata(deconv,2,'gaussian',5);

% Shift lap to start on 1, not -1 or 0.
lap = lap + (abs(min(lap))*2);

% Find laps where the landmark motor is present (on) or omitted (off)
laps_landmark_on  = find(sData.landmarks.motor.on(:,1) == 1);
on_laps  = find(ismember(lap,laps_landmark_on));

deconv     = deconv(:,on_laps);
position   = pos(on_laps);

neuron = cell(numNeuron,1);


for i = 1:numNeuron
    
    pos         = zeros(length(position),1);
    binGapsPos  = 1:N_P:max(position);
    yCenter     = 1:length(binGapsPos);
    for t = 1:length(position)
        [~,ind] = min(abs(binGapsPos-position(t)));
        pos(t) = ind;
    end
    
    % calculate Grids
    posGrid = zeros(length(pos),length(unique(pos))+1);
    for l = 1: length(pos); posGrid(l,pos(l)+1) = 1; end; posGrid(:,1) = 1;
    nBinPos = size(posGrid,2)-1;
    
    signal = deconv(i,:)';
    numFolds = 5;
    
    dt = 1/31;
    [neuron(i).testFit{n},neuron(i).trainFit{n},neuron(i).param{n}] = fit_model(posGrid,dt,signal,numFolds, beta);
    
    
end

pos_llh.neuron = neuron;
save(fullfile(savePath,'pos_llh'),'-struct', 'pos_llh');
    

function [testFit,trainFit,param_mean] = fit_model(A,dt,spiketrain,numFolds,beta)

%% this code was adapted from Hardcastle et al. 2017
% This code will section the data into 5 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 4 sections, and test the model performance on the
% remaining section. This procedure will be repeated 5 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.


%% Initialize matrices and section the data for k-fold cross-validation

[~,numCol] = size(A);
sections = numFolds*5;
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
% divide the data up into 5*num_folds pieces
edges = round(linspace(1,numel(spiketrain)+1,sections+1));

% initialize matrices
testFit = nan(numFolds,8); % var ex, correlation, llh increase, log_llh_test, log_llh_rand_test, mse, # of spikes, length of test data
trainFit = nan(numFolds,8); % var ex, correlation, llh increase, log_llh_train, log_llh_rand_train, mse, # of spikes, length of train data
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
for k = 1:numFolds
%     fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
    
    % get test data from edges - each test data chunk comes from entire session
    test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1]   ;
    
    test_spikes = spiketrain(test_ind); %test spiking
    smooth_spikes_test = conv(test_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_test = smooth_spikes_test./dt;
    test_A = A(test_ind,:);
    
    % training data
    train_ind = setdiff(1:numel(spiketrain),test_ind);
    train_spikes = spiketrain(train_ind);
    smooth_spikes_train = conv(train_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_train = smooth_spikes_train./dt;
    train_A = A(train_ind,:);
    
    opts = optimset('Gradobj','off','Hessian','off','Display','off');
    
    data{1} = train_A; data{2} = train_spikes;
    if k == 1
        init_param = 1e-3*randn(numCol, 1);
    else
        init_param = param;
    end
    
    init_param(isinf(init_param)) = 0; init_param(isnan(init_param)) = 0;
    [param] = fminunc(@(param) ln_poisson_model(param,beta,data),init_param,opts);
     
    %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_test = exp(test_A * param)/dt;
    smooth_fr_hat_test = conv(fr_hat_test,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
    sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % compute correlation
    correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson');
    
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
    r = exp(test_A * param); n = test_spikes; meanFR_test = nanmean(test_spikes); 
        
    llh_test_model      = exp(n.*log(r) - r - gammaln(n+1));%   exp(val * log(f(i,j)) - f(i,j) - gammaln(val+1))
    log_llh_test_model  = nansum(-log(llh_test_model))/sum(n);
    
    llh_test_mean       = exp(n.*log(meanFR_test) - meanFR_test - gammaln(n+1));
    log_llh_test_mean   = nansum(-log(llh_test_mean))/sum(n);
    
    log_llh_test = (-log_llh_test_model + log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test;
    
    % compute MSE
    mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
    
    % fill in all the relevant values for the test fit cases
    testFit(k,:) = [varExplain_test correlation_test log_llh_test log_llh_test_model log_llh_test_mean mse_test sum(n) numel(test_ind) ];
    
    %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_train = exp(train_A * param)/dt;
    smooth_fr_hat_train = conv(fr_hat_train,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_train-smooth_fr_train).^2);
    sst = sum((smooth_fr_train-mean(smooth_fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % compute correlation
    correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson');
    
    % compute log-likelihood
    r_train = exp(train_A * param); n_train = train_spikes; meanFR_train = nanmean(train_spikes);   
    
    llh_train_model     = exp(n_train.*log(r_train) - r_train - gammaln(n_train));
    log_llh_train_model  = nansum(-log(llh_train_model))/sum(n_train);
    
    llh_train_mean      = exp(n_train.*log(meanFR_train) - meanFR_train - gammaln(n_train));
    log_llh_train_mean  = nansum(-log(llh_train_mean))/sum(n_train);
   
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train;
    
    % compute MSE
    mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
    
    trainFit(k,:) = [varExplain_train correlation_train log_llh_train log_llh_train_model log_llh_train_mean mse_train sum(n_train) numel(train_ind)];

    % save the parameters
    paramMat(k,:) = param;

end

param_mean = nanmean(paramMat);

return
end

function [f, df, hessian] = ln_poisson_model(param,b,data)

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% compute the firing rate
u = X * param;
rate = exp(u);

% roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data

% cueOne computing the Hessian
rX = bsxfun(@times,rate,X);       
hessian_glm = rX'*X;

%% find the P, H, S, or T parameters and compute their roughness penalties

% initialize parameter-relevant variables
J = 0; J_g = []; J_h = []; 

% find the parameters
% compute the contribution for f, df, and the hessian
[J,J_g,J_h] = penalty(param,b);

%% compute f, the gradient, and the hessian 
f = sum(rate-Y(:).*u) + J ;
df = real(X' * (rate - Y(:)) + J_g);
hessian = hessian_glm + blkdiag(J_h);


end

function [J,J_g,J_h] = penalty(param,beta)
    %% smoothing functions called in the above script
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    J = beta*0.5*param'*DD1*param;
    J_g = beta*DD1*param;
    J_h = beta*DD1;
end
   
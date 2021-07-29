% Plots results 
clc,clear
rng('default');                                  % for reproducibility
addpath('results')                               % load simulation results
load('probs.mat')
load('results_bayes_net1.mat')
% --------------------- Pedro Saa UC 2021 ---------------------------------
% Extract data
thetaSMC = population_cwDK.SMCparams;
weights  = population_cwDK.weights{end};
N        = numel(weights);
indexes       = getParamIndex(probs);
numDirParams  = size(indexes,1);

% Extract true parameters from the original model
trueParamsCell = reshapeProbs(probs,probs,numDirParams,'sampling');
trueParams = [];
for ix = 1:numDirParams
    trueParams = [trueParams;(trueParamsCell{ix}(:))];
end

% Generate a random sample using the empirical sampling fxn for the SMC
% sample (get 5e4 samples for plotting)
idxSMC = [];
for ix = 1:10
    idxSMC = [idxSMC;sum(repmat(rand(N,1),1,N+1) >= cumsum([zeros(N,1),repmat(weights',N,1)],2),2)];
end

% Plot ABC marginals
plot(trueParams,1:numel(trueParams),'or','MarkerFaceColor','r')
hold on
boxplot(thetaSMC(idxSMC,:),'BoxStyle','filled','MedianStyle','line','Symbol','','Orientation','horizontal')
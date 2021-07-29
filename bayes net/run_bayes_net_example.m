% Bayesian net model
%------------------------- Pedro Saa UC 2021 ------------------------------
clc,clear

% Initialize model parameters
rng('default');                                                              % for reproducibility
load('probs.mat');

% Generate synthetic data
model.n        = 1e3;                                                        % dataset size
model.probCond = probs;                                                      % conditional probabilities
model.xdata    = simulateBayesNet(model.probCond,model.n);

% Corrupt data with blank data
noiseLevel = .3;
ix = sum(repmat(rand(numel(model.xdata),1),1,3) >= cumsum([zeros(numel(model.xdata),1),noiseLevel*ones(numel(model.xdata),1),(1-noiseLevel)*ones(numel(model.xdata),1)],2),2);
model.xmissing = model.xdata;
model.xmissing(ix==1) = 0;
model.Sdata = getSummaryStatistics(model.xmissing);

% Find full entries and use them to define the prior probabilities
fullData = model.xmissing(all(model.xmissing,2),:);

% Initialize sampling parameters 
priors    = initializePriors(fullData);
alpha     = 50;
N         = 5e3;
tolStart  = Inf;
tolFinal  = 2e-2;

%% Example 1: Use full dataset
% (1) Run example with component-wise Dirichlet kernel (MLE tuning method)
% kernel = 'cw_DK';
% method = 1;
% population_cwDK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net1','population_cwDK');
% (2) Run example with multi-component Dirichlet kernel (MLE tuning method)
% kernel = 'mc_DK';
% method = 1;
% population_mcDK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net2','population_mcDK');
% (3) Run example with multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'MGK';
% method = [];
% population_MGK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net3','population_MGK');
% (4) Run example with multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'ilrt_MGK';
% method = [];
% population_ilrtMGK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net4','population_ilrtMGK');
% (5) Run example with beta transformation and multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'beta_MGK';
% method = [];
% population_betaMGK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net5','population_betaMGK');
% (6) Run example with beta-logit transformation and multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'beta_logitMGK';
% method = [];
% population_beta_logitMGK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_bayes_net6','population_beta_logitMGK');
% (7) Run example using log-gamma and multi-variate Gaussian kernel (twice the weighted sample covariance)
kernel = 'gamma_MGK';
method = [];
population_gammaMGK = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
save('results_bayes_net7','population_gammaMGK');
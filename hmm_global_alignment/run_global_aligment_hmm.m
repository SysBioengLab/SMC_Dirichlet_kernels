% Global alignment models
%--------------------- Pedro Saa UC 2021 ----------------------------------
clc,clear

% Initialize model parameters
rng('default');                                     % for reproducibility
p0     = [.2,.5,.2,.1];                             % initial probability parameters
gamma  = .5;                                        % define model parameters           
delta  = .2;
nu     = .8;
epsil  = .1;
tau    = .1;
EMIS   = [.25,zeros(1,5);...                                                % emission probability matrix
          0,.5,.3,.15,.05,0;...
          zeros(1,5),.25;...
          zeros(1,6)];
TRANS  = [epsil,nu,0,tau;delta,gamma,delta,tau;...                          % transition probability matrix
          0,nu,epsil,tau;zeros(1,4)];

% Build probability matrices and generate synthetic data
model.n     = 1e4;                                                          % length of the sequence
model.TRANS = @(trans) [0,p0;zeros(size(trans,1),1),trans];
model.EMIS  = [zeros(1,size(EMIS,2));EMIS];
model.xdata = hist(hmmgenerate(model.n,model.TRANS(TRANS),model.EMIS));

% Initialize sampling parameters
priors{1} = ones(1,2);                                                      % prior definitions
priors{2} = ones(1,3);
alpha     = 50;
N         = 5e3;
tolStart  = Inf;
tolFinal  = 2e2;

%% Example 1: Use full dataset
% (1) Run example with component-wise Dirichlet kernel (MLE tuning method)
% kernel = 'cw_DK';
% method = 1;
% population_cwDK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full1','population_cwDK');
% (2) Run example with multi-component Dirichlet kernel (MLE tuning method)
% kernel = 'mc_DK';
% method = 1;
% population_mcDK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full2','population_mcDK');
% (3) Run example with multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'MGK';
% method = [];
% population_MGK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full3','population_MGK');
% (4) Run example with multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'ilrt_MGK';
% method = [];
% population_ilrtMGK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full4','population_ilrtMGK');
% (5) Run example with beta transformation and multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'beta_MGK';
% method = [];
% population_betaMGK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full5','population_betaMGK');
% (6) Run example with beta transformation and multi-variate Gaussian kernel (twice the weighted sample covariance)
% kernel = 'beta_logitMGK';
% method = [];
% population_beta_logitMGK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
% save('results_hmm_full6','population_beta_logitMGK');
% (7) Run example using log-gamma and multi-variate Gaussian kernel (twice the weighted sample covariance)
kernel = 'gamma_MGK';
method = [];
population_gammaMGK = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method);
save('results_hmm_full7','population_gammaMGK');
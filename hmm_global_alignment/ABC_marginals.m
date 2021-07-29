% Plots results 
clc,clear
rng('default');                                  % for reproducibility
load('results_hmm_full4')                        % load simulation results
nbins = 31;
%--------------------- Pedro Saa UC 2021 ----------------------------------
%% I) Plot ABC-posterior marginals
thetaSMC = population_ilrtMGK.HMMparams;            % pick solution from method 1
weights  = population_ilrtMGK.weights{end};

% Plot ABC marginals
figure (1)
subplot(2,3,1)
[f,x] = whistogram(thetaSMC(:,1),weights,nbins);
bar(x,f)
subplot(2,3,2)
[f,x] = whistogram(thetaSMC(:,2),weights,nbins);
bar(x,f)
subplot(2,3,3)
[f,x] = whistogram(thetaSMC(:,3),weights,nbins);
bar(x,f)
subplot(2,3,4)
[f,x] = whistogram(thetaSMC(:,4),weights,nbins);
bar(x,f)
subplot(2,3,5)
[f,x] = whistogram(thetaSMC(:,5),weights,nbins);
bar(x,f)
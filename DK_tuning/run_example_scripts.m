% Driver for the DK-tuning study behind Fig. 3 (settings in Supplementary Table S2).
%
% Runs a full factorial sweep: 4 multinomial models (3, 4, 5 and 6 components,
% one per settings<k>.mat) x 3 datasets x 3 tuning methods, i.e. 36 SMC runs.
% Within a model the dataset is enlarged in place by scaling the observed counts
% 1-, 5- and 10-fold, and the tolerance tracks it at 10% of the total count.
% The third argument of multinomial_example selects the tuning method:
%   1 MLE concentration, 2 minimum KL divergence, 3 total variance.
%
% Writes tunning_results_<k>.mat, one per model. Three conventions there are
% load-bearing and must not be "tidied":
%   - the populations are named population<method><dataset>, which is how
%     make_figures_dk.m finds them;
%   - the save is of the WHOLE workspace, because ABC_solution_comparison.m
%     also reads prior and xdata back out of tunning_results_1.mat;
%   - the file name is misspelled "tunning", and both readers expect it.
% The four blocks below are deliberately written out rather than looped, so
% that each saved file corresponds to one visible block of this script.
% --------------------- Pedro Saa UC 2026 ---------------------------------
addpath('data')

% Run 3 tuning methods with three different data sets and 3 components
clear
load settings1.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % original counts
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % counts x5; tolerance stays at 10% of the total
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % counts x10 (twice the x5 set)
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_1

% Run 3 tuning methods with three different data sets and 4 components
clear
load settings2.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % original counts
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % counts x5; tolerance stays at 10% of the total
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % counts x10 (twice the x5 set)
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_2

% Run 3 tuning methods with three different data sets and 5 components
clear
load settings3.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % original counts
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % counts x5; tolerance stays at 10% of the total
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % counts x10 (twice the x5 set)
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_3

% Run 3 tuning methods with three different data sets and 6 components
clear
load settings4.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % original counts
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % counts x5; tolerance stays at 10% of the total
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % counts x10 (twice the x5 set)
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_4
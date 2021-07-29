% Compares ABC vs analytical solution for multinomial model dataset D_33
clc,clear
rng('default');                                  % for reproducibility
addpath('results')                               % load simulation results
load('tunning_results_1.mat')
nbins = 51;
% --------------------- Pedro Saa UC 2021 ---------------------------------
%% I) Analytic posterior
% For a multinomial model, the posterior is dirichlet and the marginals are
% beta distributed
posterior = prior + xdata;
x = linspace(0,1,5e2);

% Plot analytical solution
figure (1)
subplot(1,3,1)
y1 = betapdf(x,posterior(1),sum(posterior)-posterior(1));
plot(x,y1,'-k')
xlabel('\theta_1')
ylabel('Probability')
hold on
subplot(1,3,2)
y2 = betapdf(x,posterior(2),sum(posterior)-posterior(2));
plot(x,y2,'-k')
xlabel('\theta_2')
hold on
subplot(1,3,3)
y3 = betapdf(x,posterior(3),sum(posterior)-posterior(3));
plot(x,y3,'-k')
xlabel('\theta_3')
hold on

%% II) ABC posterior approximation
thetaSMC = population13.thetaSMC;             % pick solution from method 1
weights  = population13.weights{end};

% Plot ABC marginals
figure (1)
subplot(1,3,1)
[f,x] = whistogram(thetaSMC(:,1),weights,nbins);
bar(x,f)
subplot(1,3,2)
[f,x] = whistogram(thetaSMC(:,2),weights,nbins);
bar(x,f)
subplot(1,3,3)
[f,x] = whistogram(thetaSMC(:,3),weights,nbins);
bar(x,f)
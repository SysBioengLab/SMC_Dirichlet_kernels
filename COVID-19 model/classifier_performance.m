function T = classifier_performance(resultsDir,nValidate)
% Posterior predictive performance of the COVID-19 Bayes net used as a classifier.
%
% Inputs:   resultsDir  folder holding results_bayes_net1.mat (default ./results)
%           nValidate   subjects cross-checked against BNT (default 25, 0 skips)
% Output:   T           table of indicators with 95% credible intervals
%
% For every particle of the DK posterior the predictive probability of a positive
% RT-qPCR is evaluated in closed form per subject (see predictCovidProb).
%
% The operating threshold is fixed ONCE, on the posterior-mean probability, at the
% point maximising Youden's J = sensitivity + specificity - 1, since a deployed
% classifier has a single threshold. It is then applied to every particle, so the
% reported intervals reflect posterior parameter uncertainty: point estimates are
% posterior means and the interval is the weighted 2.5th-97.5th percentile range.
% The AUC is threshold-free. Sensitivity and recall are the same quantity and are
% reported once.
%
% Table 1 of the paper is T (Youden threshold), also written to
% classifier_covid.csv. Supplementary Table S6 adds T50 (fixed 0.5 threshold),
% which is saved in classifier_covid.mat alongside it.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1) || isempty(resultsDir); resultsDir = './results'; end
if (nargin<2); nValidate = 25; end

%% Data and posterior
d     = load('xdata.mat');
xdata = d.xdata + 1;                                                       % convention false==1, true==2
yTrue = (xdata(:,1)==2);                                                   % positive class = state 2
sym   = xdata(:,2:8);
odo   = xdata(:,9:14);
n     = size(xdata,1);

f = fullfile(resultsDir,'results_bayes_net1.mat');
if ~isfile(f)
    error('classifier_performance:noResults','%s not found. Run run_bayes_net(1) first.',f);
end
S     = load(f);
pop   = S.population_cwDK;
theta = pop.SMCparams;
w     = pop.weights{end};  w = w(:)/sum(w);
N     = size(theta,1);
ESS   = 1/sum(w.^2);

off = getCovidOffsets(initializePriors_final_sim(xdata));

fprintf('DK posterior: %d particles, ESS = %.0f, n = %d subjects (%.1f%% positive)\n',...
        N,ESS,n,100*mean(yTrue));

%% Predictive probability of every subject, under every particle
P = zeros(N,n);
for i = 1:N
    P(i,:) = predictCovidProb(theta(i,:),off,sym,odo);
end

if (nValidate>0)
    validateAgainstBNT(theta(1,:),off,xdata,P(1,:),nValidate);
end

%% Operating threshold, from the posterior-mean predictive probability
pBar = (w'*P)';
[tau,J] = youdenThreshold(pBar,yTrue);
fprintf('Youden-optimal threshold: %.4f (J = %.4f)\n',tau,J);

%% Indicators, per particle, at the selected threshold and at the 0.5 default
nPos = sum(yTrue); nNeg = n-nPos;
auc  = zeros(N,1);
for i = 1:N
    auc(i) = aucRank(P(i,:)',yTrue,nPos,nNeg);
end
[T,V]   = indicatorsAt(P,yTrue,w,tau,auc,nPos,nNeg);
T50     = indicatorsAt(P,yTrue,w,0.5,auc,nPos,nNeg);

fprintf('\nOperating point: Youden-optimal threshold = %.4f\n',tau);
disp(T);
fprintf('Reference: fixed 0.5 threshold\n');
disp(T50);

acc = V(:,1); sens = V(:,2); spec = V(:,3); prec = V(:,4);
save(fullfile(resultsDir,'classifier_covid.mat'),'T','T50','tau','J','pBar',...
     'acc','sens','spec','prec','auc','w');
writetable(T,fullfile(resultsDir,'classifier_covid.csv'));
fprintf('Wrote %s and %s\n',fullfile(resultsDir,'classifier_covid.mat'),...
                            fullfile(resultsDir,'classifier_covid.csv'));

%% ROC of the posterior-mean predictive probability
[fpr,tpr] = rocCurve(pBar,yTrue);
fg = figure('Name','COVID classifier ROC','Position',[100 100 560 520],'Color','w');
plot(fpr,tpr,'-','LineWidth',1.8,'Color',[.20 .40 .70]); hold on
plot([0 1],[0 1],'k:','LineWidth',0.75)
yh = (pBar>=tau);
plot(sum(yh & ~yTrue)/nNeg,sum(yh & yTrue)/nPos,'o','MarkerSize',8,...
     'MarkerFaceColor',[.85 .33 .27],'MarkerEdgeColor','none');
hold off; axis square; grid on
xlabel('1 - specificity'); ylabel('sensitivity')
title(sprintf('ROC of the posterior-mean probabilities (AUC = %.3f); marker = Youden point',...
              aucRank(pBar,yTrue,nPos,nNeg)))                                % the AUC of THIS curve, not the posterior mean of the per-particle AUCs
legend({'ROC','chance','operating point'},'Location','southeast','Box','off')
savefig(fg,fullfile(resultsDir,'classifier_covid_roc.fig'));
print(fg,fullfile(resultsDir,'classifier_covid_roc.png'),'-dpng','-r200');
end

%% ------------------------------------------------------------------------
function [T,V] = indicatorsAt(P,yTrue,w,tau,auc,nPos,nNeg)
% Posterior mean and weighted 95% credible interval of every indicator, at a
% fixed threshold. The AUC does not depend on the threshold and is passed in.
N = size(P,1);
n = nPos + nNeg;
V = zeros(N,4);
for i = 1:N
    yh = (P(i,:)' >= tau);
    tp = sum(yh & yTrue);
    fp = sum(yh & ~yTrue);
    V(i,:) = [(tp + (nNeg-fp))/n, tp/nPos, (nNeg-fp)/nNeg, tp/max(tp+fp,1)];
end
vals  = [V,auc];
names = {'Accuracy';'Sensitivity (recall)';'Specificity';'Precision (PPV)';'AUC'};
est = zeros(5,1); lo = zeros(5,1); hi = zeros(5,1);
for j = 1:5
    est(j) = w'*vals(:,j);
    q      = wQuantile(vals(:,j),w,[0.025 0.975]);
    lo(j)  = q(1);  hi(j) = q(2);
end
T = table(names,100*est,100*lo,100*hi,'VariableNames',...
          {'indicator','estimate','ci_low','ci_high'});
end

%% ------------------------------------------------------------------------
function validateAgainstBNT(th,~,xdata,pFast,nSub)
% Cross-check the closed-form predictions against BNT's junction-tree engine on
% the first nSub subjects, using the same parameter vector. BNT needs the CPDs
% as tables rather than as the flat vector used by predictCovidProb, so the
% offsets of the caller are not needed here and the argument is ignored.
try
    addpath('./bnt-master/')
    addpath(genpathKPM('./bnt-master/'))
    L = load('Bayes_net.mat');
    bnet = L.bnet;

    priors    = initializePriors_final_sim(xdata);
    indexes   = getParamIndex(priors);
    dirBlocks = cell(size(indexes,1),1);
    for k = 1:size(indexes,1)
        dirBlocks{k} = th(indexes(k,1):indexes(k,2));
    end
    cpd = reshapeProbs(dirBlocks,priors,numel(priors),'simulation');
    for k = 1:numel(priors)
        bnet.CPD{k} = tabular_CPD(bnet,k,cpd{k}(:)');
    end
    engine = jtree_inf_engine(bnet);

    nSub = min(nSub,size(xdata,1));
    pRef = zeros(1,nSub);
    for t = 1:nSub
        ev        = cell(1,15);
        ev(2:8)   = num2cell(xdata(t,2:8));                                % symptoms
        ev(10:15) = num2cell(xdata(t,9:14));                               % odours; node 9 (sKor) stays hidden
        engine    = enter_evidence(engine,ev);
        marg      = marginal_nodes(engine,1);
        pRef(t)   = marg.T(2);
    end
    err = max(abs(pRef - pFast(1:nSub)));
    if (err < 1e-10)
        fprintf('BNT cross-check on %d subjects: PASS (max |difference| = %.2e)\n',nSub,err);
    else
        fprintf(2,'BNT cross-check on %d subjects: FAIL (max |difference| = %.2e)\n',nSub,err);
    end
catch err
    fprintf(2,'BNT cross-check skipped (%s): %s\n',err.identifier,err.message);
end
end

%% ------------------------------------------------------------------------
function a = aucRank(s,y,nPos,nNeg)
% AUC through the Mann-Whitney statistic, with mid-ranks for ties.
r = tiedrank(s);
a = (sum(r(y)) - nPos*(nPos+1)/2)/(nPos*nNeg);
end

%% ------------------------------------------------------------------------
function [fpr,tpr] = rocCurve(p,y)
[s,ix] = sort(p(:),'descend');
yy     = y(ix);
tpr    = cumsum(yy)/sum(yy);
fpr    = cumsum(~yy)/sum(~yy);
keep   = [diff(s)~=0; true];                                               % one point per distinct threshold
fpr    = [0; fpr(keep)];
tpr    = [0; tpr(keep)];
end

%% ------------------------------------------------------------------------
function v = wQuantile(x,w,q)
% Quantiles of a weighted sample, on the cumulative-midpoint convention.
[x,ix] = sort(x(:));
w      = w(ix); w = w/sum(w);
c      = cumsum(w) - 0.5*w;
keep   = [true; diff(c)>0];
v      = interp1(c(keep),x(keep),q,'linear','extrap');
v      = min(max(v,x(1)),x(end));
end

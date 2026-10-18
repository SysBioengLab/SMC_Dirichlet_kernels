function make_figure_S1(outFile)
% Exact ABC-posterior against the analytical posterior (Supplementary Fig. S1).
%
% Input:  outFile  destination .png (default ./results/FigS1_exact_posterior.png)
%
% The counts are sufficient, so running to a tolerance of exactly zero accepts
% only simulations that reproduce the observed counts: the ABC posterior then
% coincides with Dir(prior + ydata) instead of approximating it, and its marginals
% are Beta(alpha_i, alpha_0 - alpha_i).
%
% This is a separate short run, because the tuning study of Fig. 3 stops at 10% of
% the total count rather than at zero.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1) || isempty(outFile); outFile = './results/FigS1_exact_posterior.png'; end
addpath('data')

S      = load('settings1.mat');                                            % 3-component model
prior  = S.prior;
N      = S.N;
alpha  = S.alpha;
sets   = {S.xdata, 5*S.xdata};                                             % datasets D11 and D12
tags   = {'D_{11}','D_{12}'};
nbins  = 41;
D      = numel(prior);

f = figure('Name','Exact ABC posterior','Position',[60 60 1200 640],'Color','w');
for d = 1:2
    xdata = sets{d};
    rng('default');
    pop   = multinomial_example(xdata,prior,N,alpha,S.tolStart,0,1);       % method 1 (MLE), tolerance zero
    w     = pop.weights{end}; w = w(:)/sum(w);
    post  = prior + xdata;
    a0    = sum(post);
    fprintf('%s: final tolerance %g, ESS %.0f\n',tags{d},pop.tolSMC(end),1/sum(w.^2));

    for i = 1:D
        subplot(2,D,(d-1)*D+i)
        hi = betainv(0.999,post(i),a0-post(i));
        xg = linspace(0,min(1,1.15*hi),400);
        [c,dens] = wHist(pop.thetaSMC(:,i),w,nbins);
        bar(c,dens,1,'FaceColor',[.45 .55 .75],'EdgeColor','none'); hold on
        plot(xg,betapdf(xg,post(i),a0-post(i)),'k-','LineWidth',2); hold off
        xlim([0 min(1,1.15*hi)])
        xlabel(sprintf('\\theta_%d',i)); ylabel('density')
        if (i==1)
            title(sprintf('\\bf%s\\rm  dataset %s,  y = [%s]',char('A'+d-1),tags{d},...
                  strjoin(compose('%g',xdata),', ')),'HorizontalAlignment','left')
            legend({'ABC posterior','analytical'},'Box','off','FontSize',7)
        end
    end
end

d = fileparts(outFile);
if ~isempty(d) && ~isfolder(d); mkdir(d); end
print(f,outFile,'-dpng','-r200');
savefig(f,strrep(outFile,'.png','.fig'));
fprintf('Wrote %s\n',outFile);
end

%% ------------------------------------------------------------------------
function [centres,dens] = wHist(x,w,nbins)
% Weighted histogram as a density.
x  = x(:); w = w(:)/sum(w);
lo = min(x); hi = max(x);
if (hi<=lo); hi = lo+max(eps(lo),1e-12); end
edges   = linspace(lo,hi,nbins+1);
idx     = discretize(x,edges);
idx(isnan(idx)) = nbins;
dens    = accumarray(idx,w,[nbins 1]);
width   = edges(2)-edges(1);
dens    = dens/width;
centres = edges(1:end-1)+width/2;
end

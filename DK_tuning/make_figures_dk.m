function T = make_figures_dk(resultsDir)
% Post-processing for the multinomial DK-tuning study (Fig. 3).
%
% Reads tunning_results_{1..4}.mat from resultsDir (default ./results). Each
% file holds twelve runs named population<method><dataset>, with
%   method  1 MLE concentration, 2 min. KL divergence, 3 total variance
%   dataset 1 original counts, 2 five-fold, 3 ten-fold
% and corresponds to a model with 3, 4, 5 or 6 components respectively.
%
% Writes metrics_dk.mat / metrics_dk.csv and draws, each as .fig, .png and .svg,
%   fig3A_performance   cputime/ESS per model, dataset and tuning method
%   fig3B_concentration trajectory of the DK concentration parameter
%
% Fig. 3C (ABC vs analytical marginals) is produced by ABC_solution_comparison.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1); resultsDir = './results'; end

methodNames = {'MLE','min. KL','total var.'};
nComponents = [3 4 5 6];

M = [];
K = {};                                                                    % concentration trajectories
for f = 1:4
    fname = fullfile(resultsDir,sprintf('tunning_results_%d.mat',f));
    if ~isfile(fname); continue; end
    S = load(fname);
    for method = 1:3
        for dataset = 1:3
            vname = sprintf('population%d%d',method,dataset);
            if ~isfield(S,vname); continue; end
            pop = S.(vname);
            m   = getPopulationMetrics(pop,methodNames{method});
            m.model      = nComponents(f);
            m.methodIdx  = method;
            m.dataset    = dataset;
            M = [M;m];                                                     %#ok<AGROW>
            K{end+1} = struct('model',nComponents(f),'method',method,...    %#ok<AGROW>
                              'dataset',dataset,'K',pop.Kparams(:)');
        end
    end
end

if isempty(M)
    error('make_figures_dk:noResults',...
          'No tunning_results_*.mat found in %s.',resultsDir);
end

%% Metrics table
T = struct2table(M);
T = movevars(T,{'model','dataset','methodIdx','label'},'Before',1);
disp(T);
save(fullfile(resultsDir,'metrics_dk.mat'),'T','M');
writetable(T,fullfile(resultsDir,'metrics_dk.csv'));
fprintf('Wrote %s and %s\n',fullfile(resultsDir,'metrics_dk.mat'),...
                            fullfile(resultsDir,'metrics_dk.csv'));

%% Fig. 3A -- computational performance per model
models = unique(T.model)';
f1 = figure('Name','Fig 3A - performance','Position',[100 100 1200 700]);
for i = 1:numel(models)
    subplot(2,2,i)
    P = nan(3,3);                                                          % dataset x method
    for d = 1:3
        for mth = 1:3
            hit = (T.model==models(i)) & (T.dataset==d) & (T.methodIdx==mth);
            if any(hit); P(d,mth) = T.perf(hit); end
        end
    end
    bar(P)
    set(gca,'XTickLabel',{'D_{i1}','D_{i2}','D_{i3}'},'YScale','log')
    ylabel('cputime / ESS')
    title(sprintf('%d components',models(i)))
    if (i==1); legend(methodNames,'Location','northwest','Box','off'); end
    grid on
end
sgtitle('DK tuning methods: time per effective sample (lower is better)')
savefig(f1,fullfile(resultsDir,'fig3A_performance.fig'));
print(f1,fullfile(resultsDir,'fig3A_performance.png'),'-dpng','-r200');
print(f1,fullfile(resultsDir,'fig3A_performance.svg'),'-dsvg','-painters');  % vector copy for typesetting

%% Fig. 3B -- concentration trajectories for the largest dataset of each model
f2 = figure('Name','Fig 3B - concentration','Position',[100 100 1200 700]);
for i = 1:numel(models)
    subplot(2,2,i)
    hold on
    for mth = 1:3
        for k = 1:numel(K)
            if (K{k}.model==models(i)) && (K{k}.method==mth) && (K{k}.dataset==3)
                kappa = K{k}.K;
                plot(1:numel(kappa),kappa,'-o','LineWidth',1.2,'MarkerSize',3)
            end
        end
    end
    hold off
    set(gca,'YScale','log')
    xlabel('SMC iteration'); ylabel('\kappa_t')
    title(sprintf('%d components',models(i)))
    if (i==1); legend(methodNames,'Location','northwest','Box','off'); end
    grid on
end
sgtitle('DK concentration parameter along the tolerance schedule (dataset 3)')
savefig(f2,fullfile(resultsDir,'fig3B_concentration.fig'));
print(f2,fullfile(resultsDir,'fig3B_concentration.png'),'-dpng','-r200');
print(f2,fullfile(resultsDir,'fig3B_concentration.svg'),'-dsvg','-painters'); % vector copy for typesetting
fprintf('Wrote the figures into %s\n',resultsDir);
end

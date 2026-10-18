function T = getPosteriorStructure(modelName,outDir)
% Posterior dependence structure of the SMC-ABC populations, in log-ratio
% coordinates, for the kernel-comparison analysis.
%
% Inputs:   modelName  'covid', 'hmm' or 'all'
%           outDir     output folder (default ./analysis, created if missing)
% Output:   T          table of structural metrics, one row per kernel
%
% Metrics, all computed on the FINAL weighted population
%
%   Dcross    cross-block dependence index, 0.5*log( det(B) / det(Sigma) ), with
%             Sigma the posterior covariance in ILR coordinates and B its
%             block-diagonal part. This is the KL divergence from the full
%             Gaussian to its block-diagonal approximation (the trace term
%             cancels), i.e. the dependence a PRODUCT kernel cannot represent. It
%             is zero iff the blocks are independent.
%
%   condMax   largest and median within-block condition number, over blocks whose
%   condMed   ILR dimension is at least 2. The Dirichlet kernel pins its
%             covariance SHAPE to the block centre with a single free scalar, so
%             it cannot match a strongly anisotropic within-block posterior. A
%             binary node gives a 1-D block, for which anisotropy is undefined:
%             in models built mostly from binary nodes this column is vacuous and
%             the dependence story is entirely cross-block.
%
%   totVar    Aitchison total variance, sum(T_ij)/(2*d), from the variation
%             matrix T_ij = Var(log(theta_i/theta_j)). Constraint-free measure of
%             spread, and the quantity the DK's total-variance tuning is built on.
%
%   fracInfeasible / fracTolRejected / fracAccepted
%             the rejection budget over the SMC iterations, separating proposals
%             that left the simplex from those that were merely too far.
%
% Raw component correlations are deliberately NOT reported: Dirichlet variates are
% negatively correlated by the sum constraint alone, with a correlation that is a
% deterministic function of the mean. Everything here is in ILR coordinates or
% log-ratio variances, both free of that artefact.
%
% CAVEAT (circularity): each row is estimated from the population that kernel
% produced, so a poor kernel gives a poor estimate. The row flagged isReference
% (highest ESS) is the one used for the figures.
%
% A diagnostic: its tables and figures (written to outDir) back the discussion
% of posterior dependence but are not themselves items of the paper.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1) || isempty(modelName); modelName = 'all'; end
here = fileparts(mfilename('fullpath'));
if (nargin<2) || isempty(outDir); outDir = fullfile(here,'analysis'); end
if ~isfolder(outDir); mkdir(outDir); end

switch lower(modelName)
    case 'all'
        Tc = getPosteriorStructure('covid',outDir);
        Th = getPosteriorStructure('hmm',outDir);
        Tc.model = repmat({'covid'},height(Tc),1);
        Th.model = repmat({'hmm'},height(Th),1);
        T = [movevars(Tc,'model','Before',1); movevars(Th,'model','Before',1)];
        writetable(T,fullfile(outDir,'posterior_structure_all.csv'));
        compareVariationMatrices(outDir);
        return
    case 'covid'
        modelDir = fullfile(here,'COVID-19 model');
        label    = 'COVID-19 Bayes net';
    case 'hmm'
        modelDir = fullfile(here,'hmm_global_alignment');
        label    = 'pair-HMM';
    otherwise
        error('getPosteriorStructure:unknownModel',...
              'Unknown model ''%s''. Use ''covid'', ''hmm'' or ''all''.',modelName);
end

oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(modelDir);                                                         % not genpath: the COVID folder vendors the whole BNT

blockSizes = getBlockSizes(lower(modelName),modelDir);
specs      = getKernelSpec();
resultsDir = fullfile(modelDir,'results');

M    = [];
pops = {};
for k = 1:numel(specs)
    f = fullfile(resultsDir,[specs(k).resultFile,'.mat']);
    if ~isfile(f); continue; end
    S   = load(f);
    pop = S.(specs(k).varName);
    pops{end+1} = pop;                                                     %#ok<AGROW>

    m = structureMetrics(pop,blockSizes);
    m.label     = specs(k).paperName;
    m.kernelIdx = specs(k).idx;
    m.ESS       = pop.ESS(end);
    m.perf      = pop.time/pop.ESS(end);
    if isfield(pop,'rejectBudget') && (size(pop.rejectBudget,1)>1)
        b = sum(pop.rejectBudget(2:end,:),1);                              % SMC iterations only
        m.fracInfeasible   = b(2)/b(1);
        m.fracTolRejected  = b(3)/b(1);
        m.fracAccepted     = b(4)/b(1);
    else
        m.fracInfeasible = NaN; m.fracTolRejected = NaN; m.fracAccepted = NaN;
    end
    M = [M;m];                                                             %#ok<AGROW>
end

if isempty(M)
    error('getPosteriorStructure:noResults','No results found in %s.',resultsDir);
end

% Runs from different designs target different posteriors and must not be
% compared. nStats stamps the design: an unstamped run predates the current one
% and is dropped, while stamped runs that disagree are a hard error.
stamped = cellfun(@(p) isfield(p,'nStats'),pops);
if any(stamped)
    nS = cellfun(@(p) p.nStats,pops(stamped));
    if numel(unique(nS))>1
        error('getPosteriorStructure:mixedDesigns',...
              'Results use different summary sizes (%s); they are not comparable.',mat2str(unique(nS)));
    end
    if ~all(stamped)
        fprintf(2,'DROPPING %d run(s) produced before the current design: %s\n',...
                sum(~stamped),strjoin({M(~stamped).label},', '));
        M    = M(stamped);
        pops = pops(stamped);                                              %#ok<NASGU>
    end
end

[~,ref] = max([M.ESS]);
for k = 1:numel(M); M(k).isReference = (k==ref); end

T = struct2table(rmfield(M,{'variation','blockCond'}));
T = movevars(T,{'kernelIdx','label'},'Before',1);
disp(T);
save(fullfile(outDir,sprintf('posterior_structure_%s.mat',lower(modelName))),'T','M','blockSizes');
writetable(T,fullfile(outDir,sprintf('posterior_structure_%s.csv',lower(modelName))));

if strcmpi(modelName,'covid')
    fprintf(2,['NOTE: the ternary sKor node makes the COVID posterior 6-modal by label\n'...
               '      switching, and blocks 16-35 are weakly identified. Both inflate the\n'...
               '      apparent dependence; interpret Dcross for blocks 1-15 with more\n'...
               '      confidence than the aggregate.\n']);
end

makeFigures(M,ref,blockSizes,label,lower(modelName),outDir);
end

%% ------------------------------------------------------------------------
function m = structureMetrics(pop,blockSizes)
theta = pop.SMCparams;
w     = pop.weights{end};  w = w(:)/sum(w);

nClamped = sum(theta(:)<=0);
theta    = max(theta,realmin);                                             % log-ratios need strict positivity

% ILR coordinates: block-diagonal basis, one block per Dirichlet vector
Ufull = [];
for b = 1:numel(blockSizes)
    Ufull = blkdiag(Ufull,ilrBasis(blockSizes(b)));
end
Z = log(theta)*Ufull;

Sigma = weightedCov(Z,w);
Sigma = (Sigma+Sigma')/2;

% Block partition of the ILR space
dIlr   = blockSizes(:)'-1;
dIlr   = dIlr(dIlr>0);
stopI  = cumsum(dIlr);  startI = stopI-dIlr+1;

B = zeros(size(Sigma));
blockCond = nan(numel(dIlr),1);
for b = 1:numel(dIlr)
    r = startI(b):stopI(b);
    B(r,r) = Sigma(r,r);
    if (numel(r)>=2)
        ev = eig((Sigma(r,r)+Sigma(r,r)')/2);
        blockCond(b) = max(ev)/max(min(ev),eps);
    end
end

m.workDim  = size(Sigma,1);
m.nBlocks  = numel(dIlr);
m.Dcross      = 0.5*(logDet(B)-logDet(Sigma));
m.DcrossPerDim = m.Dcross/size(Sigma,1);                                   % Dcross is extensive; normalise to compare models of different size
m.condMax  = max(blockCond);
m.condMed  = median(blockCond(~isnan(blockCond)));
m.nAniso   = sum(~isnan(blockCond));                                       % blocks where anisotropy is even defined

% Aitchison variation matrix from the weighted covariance of log(theta)
C = weightedCov(log(theta),w);
d = size(C,1);
T = diag(C)*ones(1,d) + ones(d,1)*diag(C)' - 2*C;
T = max((T+T')/2,0);
T(1:d+1:end) = 0;
m.totVar    = sum(T(:))/(2*d);                                             % Aitchison total variance
m.variation = T;
m.blockCond = blockCond;
m.nClamped  = nClamped;
end

%% ------------------------------------------------------------------------
function S = weightedCov(X,w)
mu = w'*X;
dX = X - mu(ones(size(X,1),1),:);
S  = (w(:,ones(1,size(X,2))).*dX)'*dX/(1-sum(w.^2));
end

%% ------------------------------------------------------------------------
function ld = logDet(A)
% Log-determinant through the eigenvalues, with a relative ridge if the matrix
% is numerically singular (which a concentrated ABC population easily is).
ev = eig((A+A')/2);
if any(ev<=0)
    ev = ev + max(abs(ev))*1e-12;
    ev = max(ev,realmin);
end
ld = sum(log(ev));
end

%% ------------------------------------------------------------------------
function U = ilrBasis(D)
U = zeros(D,D-1);
for j = 1:D-1
    U(:,j) = sqrt(j/(j+1))*[(1/j)*ones(j,1);-1;zeros(D-j-1,1)];
end
end

%% ------------------------------------------------------------------------
function blockSizes = getBlockSizes(modelName,modelDir)
switch modelName
    case 'covid'
        d   = load(fullfile(modelDir,'xdata.mat'));
        pr  = initializePriors_final_sim(d.xdata+1);
        idx = getParamIndex(pr);
        blockSizes = (idx(:,2)-idx(:,1)+1)';
    case 'hmm'
        blockSizes = [2 3];                                                % [theta1 (2), theta2 (3)]
end
end

%% ------------------------------------------------------------------------
function makeFigures(M,ref,blockSizes,label,tag,outDir)
labels = {M.label};

% (a) variation matrix of the reference kernel
f1 = figure('Name',[label,' - variation matrix'],'Position',[100 100 720 620]);
T  = M(ref).variation;
im = imagesc(T); axis square
set(im,'AlphaData',~eye(size(T,1)));                                       % blank the zero diagonal
colorbar
title(sprintf('%s: variation matrix  T_{ij} = Var(log \\theta_i/\\theta_j)   [%s]',label,M(ref).label))
xlabel('parameter'); ylabel('parameter')
hold on
bnd = cumsum(blockSizes(:))';
for b = bnd(1:end-1)
    plot([0.5 size(T,1)+0.5],[b+0.5 b+0.5],'w-','LineWidth',0.5);
    plot([b+0.5 b+0.5],[0.5 size(T,1)+0.5],'w-','LineWidth',0.5);
end
hold off
savefig(f1,fullfile(outDir,sprintf('variation_%s.fig',tag)));
print(f1,fullfile(outDir,sprintf('variation_%s.png',tag)),'-dpng','-r200');

% (b) rejection budget
if ~all(isnan([M.fracAccepted]))
    f2 = figure('Name',[label,' - rejection budget'],'Position',[100 100 720 400]);
    bar([[M.fracInfeasible]',[M.fracTolRejected]',[M.fracAccepted]'],'stacked')
    set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',30)
    ylabel('share of proposals'); ylim([0 1])
    legend({'outside the simplex','rejected by tolerance','accepted'},'Location','eastoutside','Box','off')
    title([label,': where the proposals go'])
    savefig(f2,fullfile(outDir,sprintf('rejectbudget_%s.fig',tag)));
    print(f2,fullfile(outDir,sprintf('rejectbudget_%s.png',tag)),'-dpng','-r200');
end

% (c) cross-block dependence and within-block anisotropy
f3 = figure('Name',[label,' - dependence'],'Position',[100 100 900 380]);
subplot(1,2,1)
bar([M.Dcross],'FaceColor',[.35 .45 .70])
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',30)
ylabel('D_{cross}  (nats)'); title('Cross-block dependence a product kernel cannot represent'); grid on
subplot(1,2,2)
bc = M(ref).blockCond; bc = bc(~isnan(bc));
if isempty(bc)
    text(0.5,0.5,{'within-block anisotropy undefined:','every block is 1-D in ILR coordinates'},...
         'HorizontalAlignment','center'); axis off
else
    bar(bc,'FaceColor',[.80 .55 .30]); set(gca,'YScale','log')
    xlabel('block'); ylabel('condition number'); grid on
end
title(sprintf('Within-block anisotropy [%s]',M(ref).label))
savefig(f3,fullfile(outDir,sprintf('dependence_%s.fig',tag)));
print(f3,fullfile(outDir,sprintf('dependence_%s.png',tag)),'-dpng','-r200');
fprintf('Wrote the figures into %s\n',outDir);
end

%% ------------------------------------------------------------------------
function compareVariationMatrices(outDir)
% Side-by-side variation matrices of the two models, each from its highest-ESS
% kernel, for the visual comparison in the discussion.
f = figure('Name','Variation matrices','Position',[100 100 1250 560]);
names = {'covid','hmm'};
labs  = {'COVID-19 Bayes net','pair-HMM'};
for i = 1:2
    S = load(fullfile(outDir,sprintf('posterior_structure_%s.mat',names{i})));
    [~,ref] = max([S.M.ESS]);
    subplot(1,2,i)
    T  = S.M(ref).variation;
    im = imagesc(T); axis square
    set(im,'AlphaData',~eye(size(T,1)));
    colorbar
    dT = size(T,1);
    if (dT<=12); set(gca,'XTick',1:dT,'YTick',1:dT); end                   % integer parameter indices on small models
    xlabel('parameter'); ylabel('parameter')
    hold on
    bnd = cumsum(S.blockSizes(:))';
    for b = bnd(1:end-1)
        plot([0.5 size(T,1)+0.5],[b+0.5 b+0.5],'w-','LineWidth',0.5);
        plot([b+0.5 b+0.5],[0.5 size(T,1)+0.5],'w-','LineWidth',0.5);
    end
    hold off
    title(sprintf('%s  [%s]  D_{cross} = %.2f',labs{i},S.M(ref).label,S.M(ref).Dcross))
end
sgtitle('Variation matrices T_{ij} = Var(log \theta_i/\theta_j); white lines mark Dirichlet blocks')
savefig(f,fullfile(outDir,'variation_comparison.fig'));
print(f,fullfile(outDir,'variation_comparison.png'),'-dpng','-r200');
end

function T = make_figures_hmm(resultsDir)
% Post-processing for the pair-HMM global alignment example (Fig. 5).
%
% Reads every results file present in resultsDir (default ./results), writes the
% performance metrics to metrics_hmm.{mat,csv} and assembles Fig5.{fig,png,svg}
%   A  reserved, blank: the model schematic is added to the final artwork
%   B  cputime/ESS
%   C  effective sample size out of N
%   D  rejection distribution: outside the simplex / above tolerance / accepted
%   E  violins of the ABC-posteriors of gamma, nu and tau, one per kernel,
%      against the true values
%
% Panels B and C use a single bar colour: the kernel is named on the axis, so
% colour would be redundant. Panel E needs per-kernel colours to separate the
% five posteriors within each parameter.
%
% Palette, fonts and figure size follow figStyle, which reproduces the look of
% the final Fig. 3, and match Fig. 4 panel for panel.
%
% A run that stopped before reaching the target tolerance is included and marked
% with (*). If a kernel has no final results file but has left a checkpoint, the
% checkpoint is used and reported.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1); resultsDir = './results'; end

% True parameters of the generating model, in the order returned by the sampler:
% [gamma, delta, nu, epsilon, tau]
trueTheta = [.5,.2,.8,.1,.1];
parNames  = {'\gamma','\delta','\nu','\epsilon','\tau'};
show      = [1 3 5];                                                       % gamma, nu, tau

specs  = getKernelSpec();
M      = [];
pops   = {};
nStats = [];
for k = 1:numel(specs)
    finalFile = fullfile(resultsDir,[specs(k).resultFile,'.mat']);
    chkFile   = fullfile(resultsDir,['checkpoint_',specs(k).resultFile,'.mat']);

    if isfile(finalFile)
        S   = load(finalFile);
        pop = S.(specs(k).varName);
        tag = 'final';
    elseif isfile(chkFile)
        S   = load(chkFile);
        pop = S.population;
        tag = 'checkpoint';
    else
        continue
    end

    % Refuse to mix runs from different experimental designs. nStats stamps the
    % summary/discrepancy setup; runs predating it used a different data design
    % (one long sequence, l1 on raw counts) and are on an incomparable scale.
    if ~isfield(pop,'nStats')
        fprintf(2,'SKIPPING %s: produced before the alignment-corpus redesign; re-run this kernel.\n',...
                specs(k).resultFile);
        continue
    end
    nStats(end+1) = pop.nStats;                                            %#ok<AGROW>

    m = getPopulationMetrics(pop,specs(k).paperName);
    m.kernelIdx = specs(k).idx;
    m.kernel    = specs(k).kernel;
    m.source    = tag;
    if isfield(pop,'rejectBudget') && (size(pop.rejectBudget,1)>1)
        b = sum(pop.rejectBudget(2:end,:),1);
        m.fracInfeasible  = b(2)/b(1);
        m.fracTolRejected = b(3)/b(1);
        m.fracAccepted    = b(4)/b(1);
    else
        m.fracInfeasible = NaN; m.fracTolRejected = NaN; m.fracAccepted = NaN;
    end
    M    = [M;m];                                                          %#ok<AGROW>
    pops{end+1} = pop;                                                     %#ok<AGROW>
end

if numel(unique(nStats))>1
    error('make_figures_hmm:mixedDesigns',...
          'Results use different summary sizes (%s); they are not comparable.',mat2str(unique(nStats)));
end
if isempty(M)
    error('make_figures_hmm:noResults',...
          'No results found in %s. Run run_global_aligment_hmm(kernelIdx) first.',resultsDir);
end

%% Metrics table
T = struct2table(M);
T = movevars(T,{'kernelIdx','label','kernel'},'Before',1);
disp(T(:,{'kernelIdx','label','ESS','perf','AAR','fracInfeasible','fracAccepted','converged'}));
save(fullfile(resultsDir,'metrics_hmm.mat'),'T','M');
writetable(T,fullfile(resultsDir,'metrics_hmm.csv'));
fprintf('Wrote %s and %s\n',fullfile(resultsDir,'metrics_hmm.mat'),...
                            fullfile(resultsDir,'metrics_hmm.csv'));

labels = cellstr(T.label);                                                 % cellstr: struct2table stores a single row as char
star   = repmat({''},numel(labels),1);
star(T.converged==0) = {'*'};                                              % did not reach the target tolerance
names  = strcat(labels,star);
nK     = numel(M);

%% Assembled figure
st   = figStyle();
cols = zeros(nK,3);                                                        % one hue per kernel, panel E
for k = 1:nK
    if isKey(st.kernelColour,labels{k}); cols(k,:) = st.kernelColour(labels{k});
    else;                                cols(k,:) = st.grey;
    end
end

figW = 6.55; figH = 7.90;                                                  % inches; width matches Fig. 3
f = figure('Name','Pair-HMM','Color','w',...
           'Units','inches','Position',[0.5 0.3 figW figH]);
set(f,'PaperUnits','inches','PaperPosition',[0 0 figW figH],'PaperSize',[figW figH]);

% Panel frames, identical to Fig. 4. At this width a 12 pt title spans about
% half the page, so no row carries more than two panels, as in Fig. 3.
posA = [0.105 0.8127 0.885 0.1456];                                        % blank, reserved for the schematic
posB = [0.105 0.6076 0.375 0.1646];
posC = [0.615 0.6076 0.375 0.1646];
posD = [0.105 0.3734 0.885 0.1329];
posE = [0.105 0.0506 0.885 0.2215];

% -- A: reserved for the pair-HMM schematic --------------------------------
axA = axes('Parent',f,'Position',posA); axis(axA,'off');
panelTitle(f,posA,'\bfA',st);                                              % append the descriptive title with the schematic

% -- B: computational performance ------------------------------------------
axB = axes('Parent',f,'Position',posB);
bar(axB,T.perf,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axB,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30)
ylim(axB,[0 1.12*max(T.perf)])
yt = get(axB,'YTick'); set(axB,'YTickLabel',noLeadZero(yt))
ylabel(axB,'cputime / ESS  (s)')
styleAxes(axB,st)
panelTitle(f,posB,'\bfB\rm  Computational performance',st)

% -- C: effective sample size ----------------------------------------------
axC = axes('Parent',f,'Position',posC);
bar(axC,T.ESS,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axC,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30)
ylim(axC,[0 1.12*max(T.ESS)])
ylabel(axC,'ESS')
styleAxes(axC,st)
panelTitle(f,posC,'\bfC\rm  Effective sample size',st)

% -- D: rejection distribution ---------------------------------------------
axD = axes('Parent',f,'Position',posD);
if ~all(isnan(T.fracAccepted))
    bh = bar(axD,[T.fracInfeasible,T.fracTolRejected,T.fracAccepted],...
             'stacked','BarWidth',0.45);
    for i = 1:3; bh(i).FaceColor = st.reject(i,:); bh(i).EdgeColor = 'none'; end
    set(axD,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30,...
            'YTick',0:0.25:1,'YTickLabel',noLeadZero(0:0.25:1))
    xlim(axD,[0.4 nK+2.3])                                                 % right-hand margin holds the legend
    ylim(axD,[0 1.22])                                                     % headroom for the accepted fractions above the bars
    ylabel(axD,{'fraction of','proposals'})                                % two lines: a single one is taller than the panel
    accCol = st.reject(3,:);                                               % the accepted bar's own colour, so the two cannot drift apart
    for k = 1:nK                                                           % accepted fraction of each kernel, clear of the bar tops
        text(axD,k,1.07,sprintf('%.2f%%',100*T.fracAccepted(k)),...
             'HorizontalAlignment','center','FontName',st.font,...
             'FontSize',st.fsAnnot,'Color',accCol);
    end
    % One blue entry carries both the accepted segment and the percentages
    % printed above the bars, which are the same quantity in two renderings.
    accEntry = sprintf('\\color[rgb]{%.4f,%.4f,%.4f}accepted (%%)',accCol);
    lg = legend(axD,[bh(1) bh(2) bh(3)],...
                {'outside simplex','above tolerance',accEntry},...
                'Location','east','Box','off','FontSize',st.fsLegend,'Interpreter','tex');
    set(lg,'FontName',st.font)
else
    text(axD,0.5,0.5,'rejection budget not recorded','HorizontalAlignment','center');
    axis(axD,'off')
end
styleAxes(axD,st)
panelTitle(f,posD,'\bfD\rm  Rejection distribution',st)

% -- E: posteriors of gamma, nu and tau, one violin per kernel -------------
axE = axes('Parent',f,'Position',posE);
hold(axE,'on')
slot = 0.72/nK;                                                            % horizontal pitch between kernels
hw   = 0.40*slot;                                                          % violin half-width
hv   = gobjects(nK,1);
for g = 1:numel(show)
    p = show(g);
    for k = 1:nK
        wk = pops{k}.weights{end}; wk = wk(:)/sum(wk);
        xc = g + (k-(nK+1)/2)*slot;
        hv(k) = drawViolin(axE,xc,hw,pops{k}.HMMparams(:,p),wk,cols(k,:));
    end
    ht = plot(axE,g+[-0.44 0.44],trueTheta(p)*[1 1],'k--','LineWidth',1.1);  % true value of this parameter
    text(axE,g+0.45,trueTheta(p),sprintf(' %s = %g',parNames{p},trueTheta(p)),...
         'FontName',st.font,'FontSize',st.fsAnnot,'VerticalAlignment','middle');
end
hold(axE,'off')
set(axE,'XTick',1:numel(show),'XTickLabel',parNames(show))
xlim(axE,[0.45 numel(show)+0.72]); ylim(axE,[0.03 0.95])
yt = get(axE,'YTick'); set(axE,'YTickLabel',noLeadZero(yt))
ylabel(axE,'posterior value')
lgE = legend(axE,[hv;ht],[labels;{'true value'}],'Location','southwest',...
             'Box','off','FontSize',st.fsLegend,'NumColumns',3);            % nu sits at the top of the panel; the free space is bottom-left
set(lgE,'FontName',st.font)
styleAxes(axE,st)
panelTitle(f,posE,'\bfE\rm  ABC-posteriors of the transition probabilities',st)

savefig(f,fullfile(resultsDir,'Fig5.fig'));
print(f,fullfile(resultsDir,'Fig5.png'),'-dpng',sprintf('-r%d',st.dpi));
print(f,fullfile(resultsDir,'Fig5.svg'),'-dsvg','-painters');              % vector copy for typesetting
fprintf('Wrote the figure into %s\n',resultsDir);

%% Posterior summaries of the transition parameters, saved alongside the metrics
[~,best] = min(T.perf);                                                    % most efficient kernel available
pop      = pops{best};
w        = pop.weights{end};  w = w/sum(w);
P = table();
P.parameter = parNames(:);
P.true      = trueTheta(:);
P.postMean  = (w'*pop.HMMparams)';
P.postSD    = sqrt((w'*(pop.HMMparams-(w'*pop.HMMparams)).^2))';
disp(P);
writetable(P,fullfile(resultsDir,'posterior_hmm.csv'));
end

%% ------------------------------------------------------------------------
function panelTitle(f,pos,str,st)
% Bold letter plus title, left-aligned above the panel at the figure margin,
% as the panel titles of Fig. 3 are set.
annotation(f,'textbox',[pos(1)-0.095, pos(2)+pos(4)+0.008, pos(3)+0.095, 0.030],...
           'String',str,'EdgeColor','none','Interpreter','tex','FitBoxToText','off',...
           'Margin',0,'FontName',st.font,'FontSize',st.fsTitle,...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
end

%% ------------------------------------------------------------------------
function styleAxes(ax,st)
% Fig. 3 look: boxed axes, inward ticks, no grid, Arial throughout.
set(ax,'FontName',st.font,'FontSize',st.fsTick,'LineWidth',st.lwAxis,...
       'Box','on','TickDir','in','TickLength',[0.015 0.015],...
       'XColor','k','YColor','k','Layer','top','XGrid','off','YGrid','off');
set(get(ax,'XLabel'),'FontName',st.font,'FontSize',st.fsLabel);
set(get(ax,'YLabel'),'FontName',st.font,'FontSize',st.fsLabel);
end

%% ------------------------------------------------------------------------
function s = noLeadZero(v)
% Tick labels written without the leading zero, as in Fig. 3 (".2" not "0.2").
s = cell(numel(v),1);
for i = 1:numel(v)
    t = sprintf('%g',v(i));
    t = regexprep(t,'^(-?)0\.','$1.');
    s{i} = t;
end
end

%% ------------------------------------------------------------------------
function c = lighten(c,t)
% Blend towards white, the relation between stroke and fill colours in Fig. 3C.
c = (1-t)*c + t*[1 1 1];
end

%% ------------------------------------------------------------------------
function h = drawViolin(ax,xc,hw,x,w,col)
% One violin: a weighted kernel density mirrored about xc, with the weighted
% median marked. Width is normalised per violin, so shapes are comparable in
% form but not in area. Pale fill with a saturated outline, as in Fig. 3C.
[dens,grid] = wKde(x,w);
dens = hw*dens/max(dens);
h = patch(ax,[xc-dens; flipud(xc+dens)],[grid; flipud(grid)],lighten(col,0.45),...
          'EdgeColor',col,'LineWidth',0.9);
med = wQuantile(x,w,0.5);
plot(ax,xc,med,'o','MarkerSize',3,'MarkerFaceColor','w',...
     'MarkerEdgeColor',col,'LineWidth',0.7);
end

%% ------------------------------------------------------------------------
function [dens,grid] = wKde(x,w,nGrid)
% Gaussian kernel density of a WEIGHTED sample. The bandwidth follows
% Silverman's rule with the effective sample size in place of N, so a degenerate
% population is smoothed more, not less.
if (nargin<3); nGrid = 200; end
x  = x(:); w = w(:)/sum(w);
mu = w'*x;
sd = sqrt(max(w'*(x-mu).^2,0));
ES = 1/sum(w.^2);
h  = 1.06*sd*ES^(-1/5);
if ~(h>0); h = max(1e-9,eps(mu)); end                                      % guard a collapsed marginal
grid = linspace(min(x)-3*h,max(x)+3*h,nGrid)';
dens = exp(-0.5*((grid-x')/h).^2)*w/(h*sqrt(2*pi));
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

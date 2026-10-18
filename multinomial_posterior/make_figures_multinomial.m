function T = make_figures_multinomial(resultsDir)
% Post-processing for the multinomial benchmark with a known posterior (Fig. 4).
%
% Reads every results file present in resultsDir (default ./results) and, for
% each kernel, measures how efficiently it sampled and how closely its
% ABC-posterior reproduces the analytic posterior Dir(prior + xdata).
%
% Writes metrics_multinomial.{mat,csv} and the assembled Fig4.{fig,png,svg}
%   A  reserved, blank: the model schematic is added to the final artwork
%   B  cputime/ESS
%   C  effective sample size out of N
%   D  rejection distribution: outside the simplex / above tolerance / accepted
%   E  ABC-posterior of every component, one violin per kernel, against the
%      exact posterior
%
% Panels B and C use a single bar colour. In panel E the exact marginal
% Beta(alpha_i, alpha_0 - alpha_i) is drawn as a violin of its analytic density,
% so reference and kernels are rendered by the same convention.
%
% Palette, fonts and figure size follow figStyle, which reproduces the look of
% the final Fig. 3.
%
% The quantile summaries behind the agreement reported in the text (interquartile
% and central 95% ranges, medians, posterior means) are computed and saved even
% though the figure shows densities rather than boxes.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1); resultsDir = './results'; end

qBox = [0.025 0.25 0.50 0.75 0.975];

specs = getKernelSpec();
M     = [];
pops  = {};
for k = 1:numel(specs)
    f = fullfile(resultsDir,[specs(k).resultFile,'.mat']);
    if ~isfile(f); continue; end
    S   = load(f);
    pop = S.(specs(k).varName);
    m   = getPopulationMetrics(pop,specs(k).paperName);
    m.kernelIdx = specs(k).idx;
    m.kernel    = specs(k).kernel;
    if isfield(pop,'rejectBudget') && (size(pop.rejectBudget,1)>1)
        b = sum(pop.rejectBudget(2:end,:),1);
        m.fracInfeasible  = b(2)/b(1);
        m.fracTolRejected = b(3)/b(1);
        m.fracAccepted    = b(4)/b(1);
    else
        m.fracInfeasible = NaN; m.fracTolRejected = NaN; m.fracAccepted = NaN;
    end
    M = [M;m];                                                             %#ok<AGROW>
    pops{end+1} = pop;                                                     %#ok<AGROW>
end

if isempty(M)
    error('make_figures_multinomial:noResults',...
          'No results found in %s. Run run_multinomial(kernelIdx) first.',resultsDir);
end

postAlpha = pops{1}.posteriorAlpha;
D         = numel(postAlpha);
a0        = sum(postAlpha);
for k = 2:numel(pops)
    if ~isequal(pops{k}.posteriorAlpha,postAlpha)
        error('make_figures_multinomial:mixedPosteriors',...
              'Results were produced against different posteriors; they are not comparable.');
    end
end

% Per-kernel posterior summaries of every component, and the analytic reference
Qk    = zeros(numel(M),D,numel(qBox));                                     % kernel quantiles
Mk    = zeros(numel(M),D);                                                 % kernel posterior means
for k = 1:numel(M)
    theta = pops{k}.SMCparams;
    w     = pops{k}.weights{end};  w = w(:)/sum(w);
    for i = 1:D
        Qk(k,i,:) = wQuantile(theta(:,i),w,qBox);
        Mk(k,i)   = w'*theta(:,i);
    end
    M(k).postMean = Mk(k,:);
    M(k).meanAbsErr = mean(abs(Mk(k,:)-postAlpha/a0));                     % |E_ABC[theta] - E_exact[theta]|, averaged over components
end
Qex = zeros(D,numel(qBox));
Mex = postAlpha(:)'/a0;
for i = 1:D
    Qex(i,:) = betainv(qBox,postAlpha(i),a0-postAlpha(i));
end

%% Metrics table
T = struct2table(rmfield(M,'postMean'));
T = movevars(T,{'kernelIdx','label','kernel'},'Before',1);
disp(T(:,{'kernelIdx','label','ESS','perf','AAR','meanAbsErr','converged'}));
save(fullfile(resultsDir,'metrics_multinomial.mat'),'T','M','postAlpha','Qk','Qex','Mk','Mex');
writetable(T,fullfile(resultsDir,'metrics_multinomial.csv'));
fprintf('Wrote %s and %s\n',fullfile(resultsDir,'metrics_multinomial.mat'),...
                            fullfile(resultsDir,'metrics_multinomial.csv'));

labels = cellstr(T.label);
nK     = numel(M);

%% Assembled figure
st   = figStyle();
cols = zeros(nK,3);                                                        % one hue per kernel, panel E
for k = 1:nK
    if isKey(st.kernelColour,labels{k}); cols(k,:) = st.kernelColour(labels{k});
    else;                                cols(k,:) = st.grey;
    end
end

figW = 6.55; figH = 7.90;                                                  % inches; width matches Fig. 3 so the type prints at the same size
f = figure('Name','Multinomial benchmark','Color','w',...
           'Units','inches','Position',[0.5 0.3 figW figH]);
set(f,'PaperUnits','inches','PaperPosition',[0 0 figW figH],'PaperSize',[figW figH]);

% Panel frames. At this width a 12 pt title spans about half the page, so no
% row carries more than two panels, as in Fig. 3.
posA = [0.105 0.8127 0.885 0.1456];                                        % blank, reserved for the schematic
posB = [0.105 0.6076 0.375 0.1646];
posC = [0.615 0.6076 0.375 0.1646];
posD = [0.105 0.3734 0.885 0.1329];
posE = [0.105 0.0506 0.885 0.2215];

% -- A: reserved for the model schematic -----------------------------------
axA = axes('Parent',f,'Position',posA); axis(axA,'off');
panelTitle(f,posA,'\bfA',st);                                              % append the descriptive title with the schematic

% -- B: computational performance ------------------------------------------
axB = axes('Parent',f,'Position',posB);
bar(axB,T.perf,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axB,'XTick',1:nK,'XTickLabel',labels,'XTickLabelRotation',30,'YScale','log')
ylim(axB,[10^floor(log10(min(T.perf))), 2*max(T.perf)])                    % bars grow from the axis floor, so keep it below the smallest value
ylabel(axB,'cputime / ESS  (s)')
styleAxes(axB,st)
panelTitle(f,posB,'\bfB\rm  Computational performance',st)

% -- C: effective sample size ----------------------------------------------
axC = axes('Parent',f,'Position',posC);
bar(axC,T.ESS,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axC,'XTick',1:nK,'XTickLabel',labels,'XTickLabelRotation',30)
ylim(axC,[0 1.12*max(T.ESS)])
ylabel(axC,'ESS')
styleAxes(axC,st)
panelTitle(f,posC,'\bfC\rm  Effective sample size',st)

% -- D: rejection distribution ---------------------------------------------
axD = axes('Parent',f,'Position',posD);
if ~all(isnan([M.fracAccepted]))
    bh = bar(axD,[[M.fracInfeasible]',[M.fracTolRejected]',[M.fracAccepted]'],...
             'stacked','BarWidth',0.45);
    for i = 1:3; bh(i).FaceColor = st.reject(i,:); bh(i).EdgeColor = 'none'; end
    set(axD,'XTick',1:nK,'XTickLabel',labels,'XTickLabelRotation',30,...
            'YTick',0:0.25:1,'YTickLabel',noLeadZero(0:0.25:1))
    xlim(axD,[0.4 nK+2.3])                                                 % right-hand margin holds the legend
    ylim(axD,[0 1.22])                                                     % headroom for the accepted fractions above the bars
    ylabel(axD,{'fraction of','proposals'})                                % two lines: a single one is taller than the panel
    accCol = st.reject(3,:);                                               % the accepted bar's own colour, so the two cannot drift apart
    for k = 1:nK                                                           % accepted fraction of each kernel, clear of the bar tops
        text(axD,k,1.07,sprintf('%.2f%%',100*M(k).fracAccepted),...
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

% -- E: ABC-posterior of every component against the exact posterior -------
axE = axes('Parent',f,'Position',posE);
hold(axE,'on')
slot = 0.80/(nK+1);                                                        % pitch: the exact posterior plus one violin per kernel
hw   = 0.40*slot;
h    = gobjects(nK+1,1);
qTail = [0.005 0.995];                                                     % every violin is truncated to the same central range
yTop  = 0;
for i = 1:D
    rEx    = betainv(qTail,postAlpha(i),a0-postAlpha(i));
    gEx    = linspace(rEx(1),rEx(2),200)';
    h(1)   = drawViolin(axE,i+(1-(nK+2)/2)*slot,hw,...
                        betapdf(gEx,postAlpha(i),a0-postAlpha(i)),gEx,Qex(i,3),st.black);
    yTop   = max(yTop,rEx(2));
    for k = 1:nK
        wk      = pops{k}.weights{end}; wk = wk(:)/sum(wk);
        xs      = pops{k}.SMCparams(:,i);
        rK      = wQuantile(xs,wk,qTail);
        [dk,gk] = wKde(xs,wk,rK);
        h(k+1)  = drawViolin(axE,i+(k+1-(nK+2)/2)*slot,hw,dk,gk,Qk(k,i,3),cols(k,:));
        yTop    = max(yTop,rK(2));
    end
end
hold(axE,'off')
set(axE,'XTick',1:D,'XTickLabel',arrayfun(@(i) sprintf('\\theta_%d',i),1:D,'UniformOutput',false))
xlim(axE,[0.5 D+0.5]); ylim(axE,[0 1.16*yTop])
yt = get(axE,'YTick'); set(axE,'YTickLabel',noLeadZero(yt))
ylabel(axE,'posterior value')
lgE = legend(axE,h,[{'exact'};labels],'Location','northwest','Box','off',...
             'FontSize',st.fsLegend,'NumColumns',3);
set(lgE,'FontName',st.font)
styleAxes(axE,st)
panelTitle(f,posE,'\bfE\rm  Component-wise ABC-posterior against the exact posterior',st)

savefig(f,fullfile(resultsDir,'Fig4.fig'));
print(f,fullfile(resultsDir,'Fig4.png'),'-dpng',sprintf('-r%d',st.dpi));
print(f,fullfile(resultsDir,'Fig4.svg'),'-dsvg','-painters');              % vector copy for typesetting
fprintf('Wrote the figure into %s\n',resultsDir);
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
function h = drawViolin(ax,xc,hw,dens,grid,med,col)
% One violin: a density mirrored about xc, with the median marked. Width is
% normalised per violin, so shapes are comparable in form but not in area.
% Pale fill with a saturated outline, as the histograms of Fig. 3C.
dens = hw*dens/max(dens);
h = patch(ax,[xc-dens; flipud(xc+dens)],[grid; flipud(grid)],lighten(col,0.45),...
          'EdgeColor',col,'LineWidth',0.9);
plot(ax,xc,med,'o','MarkerSize',3,'MarkerFaceColor','w',...
     'MarkerEdgeColor',col,'LineWidth',0.7);
end

%% ------------------------------------------------------------------------
function [dens,grid] = wKde(x,w,rng,nGrid)
% Gaussian kernel density of a WEIGHTED sample, evaluated over rng = [lo hi]. The
% bandwidth follows Silverman's rule with the effective sample size in place of
% N, so a degenerate population is smoothed more, not less. The whole sample
% contributes to the density; rng only sets the window that is drawn.
if (nargin<4); nGrid = 200; end
x  = x(:); w = w(:)/sum(w);
mu = w'*x;
sd = sqrt(max(w'*(x-mu).^2,0));
ES = 1/sum(w.^2);
h  = 1.06*sd*ES^(-1/5);
if ~(h>0); h = max(1e-9,eps(mu)); end                                      % guard a collapsed marginal
if (nargin<3) || isempty(rng); rng = [min(x)-3*h, max(x)+3*h]; end
grid = linspace(max(rng(1),0),rng(2),nGrid)';                              % components are non-negative
dens = exp(-0.5*((grid-x')/h).^2)*w/(h*sqrt(2*pi));
end

%% ------------------------------------------------------------------------
function v = wQuantile(x,w,p)
% Quantiles of a weighted sample, on the cumulative-midpoint convention.
[x,ix] = sort(x(:));
w      = w(ix); w = w/sum(w);
c      = cumsum(w) - 0.5*w;
keep   = [true; diff(c)>0];                                                % interp1 needs strictly increasing nodes
v      = interp1(c(keep),x(keep),p,'linear','extrap');
v      = min(max(v,x(1)),x(end));
end

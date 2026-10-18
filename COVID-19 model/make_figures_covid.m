function T = make_figures_covid(resultsDir)
% Post-processing for the COVID-19 screening Bayes net example (Fig. 6).
%
% Reads every results file present in resultsDir (default ./results), writes the
% performance metrics to metrics_covid.{mat,csv} and assembles Fig6.{fig,png,svg}
%   A  reserved, blank: the network schematic is added to the final artwork
%   B  tolerance schedules
%   C  cputime/ESS, logarithmic
%   D  effective sample size out of N, logarithmic
%   E  rejection distribution: outside the simplex / above tolerance / accepted
%   F  ROC of the DK posterior, with the operating point at the Youden threshold
%
% Panels C and D use a single bar colour: the kernel is named on the axis, so
% colour would be redundant. Panel B needs per-kernel colours to separate five
% overlapping curves.
%
% Palette, fonts and figure size follow figStyle, which reproduces the look of
% the final Fig. 3, as in Figs. 4 and 5. This figure carries six panels rather
% than five, so the bar panels of the third row are of unequal width: the
% rejection panel has to hold its own legend as well as a label per bar.
%
% It also draws covid_marginals, a diagnostic with one weighted marginal per
% Dirichlet block; the paper has no such panel, but this is what reveals the
% label switching of the ternary sKor node.
%
% A run that stopped before reaching the target tolerance is included and marked
% with (*). If a kernel has no final results file but has left a checkpoint, the
% checkpoint is used and reported.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1); resultsDir = './results'; end

specs = getKernelSpec();
M     = [];
pops  = {};
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

if isempty(M)
    error('make_figures_covid:noResults',...
          'No results found in %s. Run run_bayes_net(kernelIdx) first.',resultsDir);
end

%% Metrics table
T = struct2table(M);
T = movevars(T,{'kernelIdx','label','kernel'},'Before',1);
disp(T(:,{'kernelIdx','label','ESS','perf','AAR','fracInfeasible','fracAccepted','converged'}));
save(fullfile(resultsDir,'metrics_covid.mat'),'T','M');
writetable(T,fullfile(resultsDir,'metrics_covid.csv'));
fprintf('Wrote %s and %s\n',fullfile(resultsDir,'metrics_covid.mat'),...
                            fullfile(resultsDir,'metrics_covid.csv'));

labels = cellstr(T.label);                                                 % cellstr: struct2table stores a single row as char
star   = repmat({''},numel(labels),1);
star(T.converged==0) = {'*'};                                              % did not reach the target tolerance
names  = strcat(labels,star);
nK     = numel(M);

data        = load('xdata.mat');
xdata       = data.xdata + 1;                                              % convention false==1, true==2
priorShapes = initializePriors_final_sim(xdata);
indexes     = getParamIndex(priorShapes);

%% Assembled figure
st   = figStyle();
cols = zeros(nK,3);                                                        % one hue per kernel, panel B
for k = 1:nK
    if isKey(st.kernelColour,labels{k}); cols(k,:) = st.kernelColour(labels{k});
    else;                                cols(k,:) = st.grey;
    end
end

% Two rows of three. Three columns cannot carry a 12 pt title at the 6.55 in
% width of Figs. 3-5, so the canvas is widened instead of the type reduced: the
% point sizes stay those of figStyle and every legend fits inside its own panel.
figW = 9.60; figH = 5.90;
f = figure('Name','COVID-19 Bayes net','Color','w',...
           'Units','inches','Position',[0.3 0.6 figW figH]);
set(f,'PaperUnits','inches','PaperPosition',[0 0 figW figH],'PaperSize',[figW figH]);

colX = [0.062 0.389 0.716];                                                % three equal columns
colW = 0.272;
rowY = [0.5966 0.1102];                                                    % top row, bottom row
rowH = 0.3390;
posA = [colX(1) rowY(1) colW rowH];                                        % blank, reserved for the schematic
posB = [colX(2) rowY(1) colW rowH];
posC = [colX(3) rowY(1) colW rowH];
posD = [colX(1) rowY(2) colW rowH];
posE = [colX(2) rowY(2) colW rowH];
posF = [colX(3) rowY(2) colW rowH];

% -- A: reserved for the network schematic ---------------------------------
axA = axes('Parent',f,'Position',posA); axis(axA,'off');
panelTitle(f,posA,'\bfA',st);                                              % append the descriptive title with the schematic

% -- B: tolerance schedules -------------------------------------------------
axB = axes('Parent',f,'Position',posB);
hold(axB,'on')
for k = 1:nK
    plot(axB,0:numel(pops{k}.tolSMC)-1,pops{k}.tolSMC,'-o','Color',cols(k,:),...
         'LineWidth',st.lwCurve,'MarkerSize',2.5,'MarkerFaceColor',cols(k,:),...
         'MarkerEdgeColor','none');
end
tolF  = min(cellfun(@(p) p.tolSMC(end),pops));
maxIt = max(cellfun(@(p) numel(p.tolSMC),pops)) - 1;
xlim(axB,[0 maxIt+0.6]); ylim(axB,ylim(axB));                              % fix both limits before the rule is drawn
xl = xlim(axB); yl = ylim(axB);
plot(axB,xl,[tolF tolF],'k--','LineWidth',0.75);                           % spans the whole panel
xlim(axB,xl);                                                              % its own endpoints must not re-expand the axis
text(axB,xl(1)+0.025*diff(xl),tolF+0.055*diff(yl),sprintf('\\epsilon_T = %g',tolF),...
     'FontName',st.font,'FontSize',st.fsAnnot);
hold(axB,'off')
xlabel(axB,'SMC iteration'); ylabel(axB,'tolerance \epsilon_t')
styleAxes(axB,st)
lgB = legend(axB,labels,'Location','northeast','Box','off','FontSize',st.fsLegend);
set(lgB,'FontName',st.font)
panelTitle(f,posB,'\bfB\rm  Tolerance schedule',st)

% -- C: computational performance -------------------------------------------
axC = axes('Parent',f,'Position',posC);
bar(axC,T.perf,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axC,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30,'YScale','log')
loC = floor(log10(min(T.perf))); hiC = ceil(log10(3*max(T.perf)));
ylim(axC,[10^loC, 10^hiC])                                                 % bars grow from the axis floor
set(axC,'YTick',10.^(loC:hiC))                                             % short axis: name every decade or MATLAB drops most
ylabel(axC,'cputime / ESS  (s)')
styleAxes(axC,st)
panelTitle(f,posC,'\bfC\rm  Computational performance',st)

% -- D: effective sample size -----------------------------------------------
axD = axes('Parent',f,'Position',posD);
bar(axD,T.ESS,'FaceColor',st.bar,'EdgeColor','none','BarWidth',0.65)
set(axD,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30,'YScale','log')
hiD = ceil(log10(max(T.ESS)));
ylim(axD,[1, 10^(hiD+0.55)])
set(axD,'YTick',10.^(0:hiD))                                               % short axis: name every decade or MATLAB drops most
for k = 1:nK                                                               % the span is three decades, so print the value
    text(axD,k,T.ESS(k)*1.7,sprintf('%.0f',T.ESS(k)),...
         'HorizontalAlignment','center','FontName',st.font,'FontSize',st.fsAnnot);
end
ylabel(axD,'ESS')
styleAxes(axD,st)
panelTitle(f,posD,'\bfD\rm  Effective sample size',st)

% -- E: rejection distribution ----------------------------------------------
axE = axes('Parent',f,'Position',posE);
if ~all(isnan(T.fracAccepted))
    bh = bar(axE,[T.fracInfeasible,T.fracTolRejected,T.fracAccepted],...
             'stacked','BarWidth',0.45);
    for i = 1:3; bh(i).FaceColor = st.reject(i,:); bh(i).EdgeColor = 'none'; end
    set(axE,'XTick',1:nK,'XTickLabel',names,'XTickLabelRotation',30,...
            'YTick',0:0.25:1,'YTickLabel',noLeadZero(0:0.25:1))
    xlim(axE,[0.4 nK+0.6])
    ylim(axE,[0 1.80])                                                     % the legend shares the panel: three columns leave no room beside it
    ylabel(axE,'fraction of proposals')
    accCol = st.reject(3,:);                                               % the accepted bar's own colour, so the two cannot drift apart
    for k = 1:nK                                                           % accepted fraction of each kernel, clear of the bar tops
        text(axE,k,1.07,sprintf('%.2f%%',100*T.fracAccepted(k)),...
             'HorizontalAlignment','center','FontName',st.font,...
             'FontSize',st.fsAnnot-1,'Color',accCol);
    end
    % One blue entry carries both the accepted segment and the percentages
    % printed above the bars, which are the same quantity in two renderings.
    accEntry = sprintf('\\color[rgb]{%.4f,%.4f,%.4f}accepted (%%)',accCol);
    lgE = legend(axE,[bh(1) bh(2) bh(3)],...
                 {'outside simplex','above tolerance',accEntry},...
                 'Location','north','Box','off',...
                 'FontSize',st.fsLegend-1,'Interpreter','tex');
    set(lgE,'FontName',st.font)                                            % a point down: this row is the tightest in the figure
else
    text(axE,0.5,0.5,'rejection budget not recorded','HorizontalAlignment','center');
    axis(axE,'off')
end
styleAxes(axE,st)
panelTitle(f,posE,'\bfE\rm  Rejection distribution',st)

% -- F: posterior ROC of the DK classifier, with its 95% credible band ------
% Every quantity here is the one reported in Table 1. The curve is the posterior
% MEAN ROC, obtained by evaluating each particle's ROC on a common grid of
% 1 - specificity and averaging vertically. The AUC quoted is the posterior mean
% of the per-particle areas, which is what Table 1 gives; the ROC of the
% posterior-mean probabilities would be a different, slightly larger number.
axF = axes('Parent',f,'Position',posF);
dk  = find(T.kernelIdx==1,1);
if isempty(dk)
    text(axF,0.5,0.5,'DK results not available','HorizontalAlignment','center');
    axis(axF,'off')
else
    yTrue = (xdata(:,1)==2);
    off   = getCovidOffsets(priorShapes);
    popDK = pops{dk};
    wDK   = popDK.weights{end}; wDK = wDK(:)/sum(wDK);
    nPart = size(popDK.SMCparams,1);
    Pdk   = zeros(nPart,size(xdata,1));
    for i = 1:nPart
        Pdk(i,:) = predictCovidProb(popDK.SMCparams(i,:),off,xdata(:,2:8),xdata(:,9:14));
    end

    nPos = sum(yTrue); nNeg = numel(yTrue)-nPos;
    fg   = linspace(0,1,201)';                                             % common 1 - specificity grid
    Tp   = zeros(nPart,numel(fg));
    auc  = zeros(nPart,1);
    for i = 1:nPart
        [fp,tp] = rocCurve(Pdk(i,:)',yTrue);
        auc(i)  = trapz(fp,tp);                                            % on this particle's own breakpoints: exact
        [fu,iu] = unique(fp,'last');                                       % interp1 needs unique nodes
        Tp(i,:) = interp1(fu,tp(iu),fg,'linear');
    end
    mTpr = (wDK'*Tp)';                                                     % posterior mean ROC
    band = zeros(numel(fg),2);
    for g = 1:numel(fg)
        band(g,:) = wQuantileCol(Tp(:,g),wDK,[0.025 0.975]);
    end
    % Integrating the vertically averaged curve instead would read ~8e-4 higher,
    % because averaging at fixed 1 - specificity is not the same operation as
    % averaging the areas; below the reported precision, but not interchangeable.
    aucM  = wDK'*auc;
    aucCI = wQuantileCol(auc,wDK,[0.025 0.975]);

    % operating point of Table 1: posterior mean sensitivity and specificity at
    % the Youden threshold, with their credible intervals as error bars
    pBar    = (wDK'*Pdk)';
    tau     = youdenThreshold(pBar,yTrue);
    yh      = (Pdk >= tau);
    sensP   = sum(yh(:,yTrue),2)/nPos;
    specP   = sum(~yh(:,~yTrue),2)/nNeg;
    sM      = wDK'*sensP;   sQ = wQuantileCol(sensP,wDK,[0.025 0.975]);
    pM      = wDK'*specP;   pQ = wQuantileCol(specP,wDK,[0.025 0.975]);

    hold(axF,'on')
    hB = fill(axF,[fg;flipud(fg)],[band(:,1);flipud(band(:,2))],...
              lighten(st.navy,0.72),'EdgeColor','none');
    hC = plot(axF,[0 1],[0 1],'k:','LineWidth',0.75);
    hM = plot(axF,fg,mTpr,'-','LineWidth',st.lwCurve,'Color',st.navy);
    hO = errorbar(axF,1-pM,sM,sM-sQ(1),sQ(2)-sM,(1-pM)-(1-pQ(2)),(1-pQ(1))-(1-pM),...
                  'o','MarkerSize',5,'MarkerFaceColor',st.magenta,...
                  'MarkerEdgeColor','none','Color',st.magenta,'LineWidth',0.9,'CapSize',3);
    hold(axF,'off')
    xlim(axF,[0 1]); ylim(axF,[0 1])
    set(axF,'XTick',0:0.2:1,'XTickLabel',noLeadZero(0:0.2:1),...
            'YTick',0:0.2:1,'YTickLabel',noLeadZero(0:0.2:1))
    xlabel(axF,'1 - specificity'); ylabel(axF,'sensitivity')
    styleAxes(axF,st)

    % Legend bottom-right, under the chance line; the AUC of Table 1 goes
    % top-left, in the corner the curve and its band leave free. Split this way
    % both fit inside the panel, which one combined entry would not.
    text(axF,0.045,0.975,{sprintf('AUC = %.1f%%',100*aucM),...
                          sprintf('[%.1f, %.1f]',100*aucCI(1),100*aucCI(2))},...
         'FontName',st.font,'FontSize',st.fsAnnot,...
         'HorizontalAlignment','left','VerticalAlignment','top');
    lgF = legend(axF,[hM hB hO hC],...
                 {'posterior mean','95% credible band',...
                  sprintf('operating point (\\tau = %.2f)',tau),'chance'},...
                 'Location','southeast','Box','off',...
                 'FontSize',st.fsLegend,'Interpreter','tex');
    set(lgF,'FontName',st.font);
    fprintf('Fig. 6F: AUC %.4f [%.4f, %.4f]; operating point tau = %.4f, sens %.4f, spec %.4f\n',...
            aucM,aucCI(1),aucCI(2),tau,sM,pM);
end
panelTitle(f,posF,'\bfF\rm  ROC of the DK classifier',st)

savefig(f,fullfile(resultsDir,'Fig6.fig'));
print(f,fullfile(resultsDir,'Fig6.png'),'-dpng',sprintf('-r%d',st.dpi));
print(f,fullfile(resultsDir,'Fig6.svg'),'-dsvg','-painters');              % vector copy for typesetting

%% Diagnostic marginals, one panel per Dirichlet block
nBlocks     = size(indexes,1);
nbins       = 31;

[~,best] = min(T.perf);                                                    % most efficient kernel available
pop      = pops{best};
w        = pop.weights{end};  w = w/sum(w);

f2 = figure('Name','COVID ABC marginals (diagnostic)','Position',[100 100 1400 800]);
nCol = 7;
nRow = ceil(nBlocks/nCol);
for b = 1:nBlocks
    subplot(nRow,nCol,b)
    [c,d] = wHist(pop.SMCparams(:,indexes(b,1)),w,nbins);
    bar(c,d,1,'FaceColor',[.45 .55 .75],'EdgeColor','none')
    xlim([0 1]); set(gca,'YTick',[])
    title(sprintf('block %d',b),'FontWeight','normal','FontSize',8)
end
sgtitle(sprintf('ABC-posterior marginals, first component of each Dirichlet block (%s)',labels{best}));

savefig(f2,fullfile(resultsDir,'covid_marginals.fig'));
print(f2,fullfile(resultsDir,'covid_marginals.png'),'-dpng','-r150');
fprintf('Wrote the figures into %s\n',resultsDir);
end

%% ------------------------------------------------------------------------
function panelTitle(f,pos,str,st)
% Bold letter plus title, left-aligned above the panel at the figure margin,
% as the panel titles of Fig. 3 are set.
annotation(f,'textbox',[pos(1)-0.055, pos(2)+pos(4)+0.010, pos(3)+0.055, 0.042],...
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
function v = wQuantileCol(x,w,q)
% Quantiles of a weighted sample, on the cumulative-midpoint convention.
[x,ix] = sort(x(:));
w      = w(ix); w = w/sum(w);
c      = cumsum(w) - 0.5*w;
keep   = [true; diff(c)>0];
if (sum(keep)<2); v = repmat(x(1),size(q)); return; end                    % degenerate column
v = interp1(c(keep),x(keep),q,'linear','extrap');
v = min(max(v,x(1)),x(end));
end

%% ------------------------------------------------------------------------
function [fpr,tpr] = rocCurve(p,y)
% ROC of a score vector: one point per distinct threshold, anchored at the origin.
[s,ix] = sort(p(:),'descend');
yy     = y(ix);
tpr    = cumsum(yy)/sum(yy);
fpr    = cumsum(~yy)/sum(~yy);
keep   = [diff(s)~=0; true];
fpr    = [0; fpr(keep)];
tpr    = [0; tpr(keep)];
end

%% ------------------------------------------------------------------------
function [centres,dens] = wHist(x,w,nbins)
% Weighted histogram as a density, written out so that a degenerate marginal
% (all particles at one value) cannot error.
x  = x(:); w = w(:)/sum(w);
lo = min(x); hi = max(x);
if (hi<=lo); hi = lo+max(eps(lo),1e-12); end
edges = linspace(lo,hi,nbins+1);
idx   = discretize(x,edges);
idx(isnan(idx)) = nbins;
dens    = accumarray(idx,w,[nbins 1]);
width   = edges(2)-edges(1);
dens    = dens/width;
centres = edges(1:end-1)+width/2;
end

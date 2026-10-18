function S = figStyle()
% Shared look of the article figures, taken from the final Fig. 3.
%
% Fig. 3 is 6.55 x 6.05 in at 300 dpi, set in Arial, with bold panel letters
% followed by a regular-weight title, no grid lines, and tick labels written
% without a leading zero. The palette below is the one sampled from that file.
%
% Colours come in pairs: a saturated stroke and the same hue lightened ~45%
% towards white for the fill, which is how the histograms of Fig. 3C are drawn.
% --------------------- Pedro Saa UC 2026 ---------------------------------

S.font      = 'Arial';
S.fsLetter  = 12;                        % bold panel letter
S.fsTitle   = 12;                        % panel title, regular weight
S.fsLabel   = 11;                        % axis labels
S.fsTick    = 10;                        % tick labels
S.fsLegend  = 9;
S.fsAnnot   = 9;                         % in-panel annotations
S.lwAxis    = 0.75;
S.lwCurve   = 1.5;                       % 1.44 pt measured on the Fig. 3C curves
S.dpi       = 300;

% -- palette sampled from Fig3.png -----------------------------------------
S.navy    = [ 46  49 146]/255;
S.green   = [ 43 182 115]/255;
S.magenta = [236   0 140]/255;
S.cyan    = [  0 174 239]/255;
S.orange  = [247 148  30]/255;
S.lime    = [140 198  63]/255;
S.grey    = [147 149 152]/255;
S.black   = [  0   0   0];

% Single bar colour for the per-kernel bar panels
S.bar     = S.navy;

% Stacked rejection categories: inadmissible / too far / accepted. The accepted
% category takes the bar-plot blue, and the accepted fractions printed above the
% bars read their colour from here, so segment and label cannot drift apart.
S.reject  = [S.magenta; S.grey; S.bar];

% One hue per kernel, keyed by the label used in the metrics table. Five
% clearly separated hues; the exact posterior is drawn in black, as the
% analytical marginal of Fig. 3C is.
S.kernelColour = containers.Map(...
    {'DK','LT-MGK','ILRT-MGK','BLT-MGK','LGT-MGK'}, ...
    {S.navy, S.cyan, S.green, S.orange, S.magenta});
end

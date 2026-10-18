function [ff,xlocs] = whistogram(xdata,ww,nbins,minX,maxX)
% Weighted histogram of xdata as a density (unit area).
%
% Inputs:   xdata  samples;  ww  their weights;  nbins  number of bin EDGES, so
%           nbins-1 bins are returned;  minX, maxX  optional range
% Outputs:  ff     (nbins-1) x 1 densities;  xlocs  bin positions
%
% xlocs uses delta = range/nbins rather than the bin width range/(nbins-1), which
% places each bar about range/(2*nbins^2) left of its true centre -- invisible at
% the bin counts used, and kept as is because Fig. 3C was drawn with it.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if nargin<4
    minX  = min(xdata);
    maxX  = max(xdata);
elseif nargin<5
    maxX  = max(xdata);
end
xlocs = linspace(minX,maxX,nbins);
delta = (maxX-minX)/nbins;

% Make sure that the weights sum to 1
ww = ww/sum(ww);

% Sort xdata
[xdata,ic] = sort(xdata,'ascend');
ww         = ww(ic);

% Compute probabilities for each bin
pw    = cumsum(ww);
ff    = zeros(numel(xlocs(:))-1,1);
for ix = 2:numel(xlocs(:))
    if ix<numel(xlocs(:))
        jx = find(xdata>xlocs(ix),1,'first');
        if (ix==2)
            ff(ix-1) = pw(jx-1);
        else
            ff(ix-1) = pw(jx-1)-pw(prevIx);
        end
    else
        ff(end) = pw(end)-pw(prevIx);
    end
    prevIx = jx-1;
end

% Center bins
xlocs = (xlocs(1:end-1)+delta/2)';
ff    = ff/trapz(xlocs,ff);

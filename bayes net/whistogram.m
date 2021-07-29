function [ff,xlocs] = whistogram(xdata,ww,nbins,minX,maxX)
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

function S = getSummaryStatistics(xdata)
% Get summary statistics
% --------------------- Pedro Saa UC 2021 ---------------------------------
m = size(xdata,2);
n = max(xdata(:));
S = zeros(m,n);
for ix = 1:m
    for jx = 1:n
        S(ix,jx) = sum(xdata(:,ix)==jx);
    end    
end
S = S./repmat(sum(S,2),1,n);
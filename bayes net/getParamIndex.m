function [indexes] = getParamIndex(priors)
% Get indexes for the parameters based on the prior
% --------------------- Pedro Saa UC 2021 ---------------------------------
indexes = [];
prevIdx = 1;
for ix = 1:numel(priors)
    for jx = 1:size(priors{ix},1)
        currIdx    = prevIdx + numel(priors{ix}(jx,:)) - 1;
        indexes    = [indexes;[prevIdx,currIdx]];
        prevIdx    = currIdx + 1;
    end
end
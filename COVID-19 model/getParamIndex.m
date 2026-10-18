function [indexes] = getParamIndex(priors)
% Column range of every Dirichlet block in the flat parameter vector.
%
% Input:    priors   cell of Dirichlet hyper-parameters; each row of priors{k} is
%                    one block (one row of the node's conditional probability table)
% Output:   indexes  nBlocks x 2 [first last] columns, blocks in node order
% --------------------- Pedro Saa UC 2026 ---------------------------------
indexes = [];
prevIdx = 1;
for ix = 1:numel(priors)
    for jx = 1:size(priors{ix},1)
        currIdx    = prevIdx + numel(priors{ix}(jx,:)) - 1;
        indexes    = [indexes;[prevIdx,currIdx]];
        prevIdx    = currIdx + 1;
    end
end
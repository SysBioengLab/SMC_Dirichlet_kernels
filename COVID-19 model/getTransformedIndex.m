function transIdxs = getTransformedIndex(indexes)
% Column ranges of each Dirichlet block in the TRANSFORMED parameter vector.
%
% Input:    indexes    nBlocks x 2 [first last] columns of each block in the
%                      simplex vector (see getParamIndex)
% Output:   transIdxs  nBlocks x 2 ranges in the (D-1)-per-block space used by the
%                      LT-, ILRT- and BLT-MGK, where each block loses one coordinate
% --------------------- Pedro Saa UC 2026 ---------------------------------
blockLen    = indexes(:,2)-indexes(:,1);
transIdxs = zeros(size(indexes));
prevIdx   = 1;
for ix = 1:numel(blockLen)
    transIdxs(ix,1) = prevIdx;
    transIdxs(ix,2) = prevIdx + blockLen(ix) - 1;
    prevIdx         = transIdxs(ix,2) + 1;
end

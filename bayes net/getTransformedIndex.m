function transIdxs = getTransformedIndex(indexes)
% Returns indexes on the transformed space
% --------------------- Pedro Saa UC 2021 ---------------------------------
length    = indexes(:,2)-indexes(:,1);
transIdxs = zeros(size(indexes));
prevIdx   = 1;
for ix = 1:numel(length)    
    transIdxs(ix,1) = prevIdx;
    transIdxs(ix,2) = prevIdx + length(ix) - 1;
    prevIdx         = transIdxs(ix,2) + 1;
end

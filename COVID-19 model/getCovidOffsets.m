function off = getCovidOffsets(priors)
% Flat-vector offsets of each node's conditional probability table.
%
% Input:    priors  15x1 cell of Dirichlet hyper-parameters, sized
%                   (parent states) x (own states)
% Output:   off     struct with fields start/rows/cols, each 1x15
%
% theta concatenates the Dirichlet blocks in node order and each node's blocks are
% consecutive, so
%   P(node k = s | parent j) = theta( off.start(k) + (j-1)*off.cols(k) + (s-1) )
% letting the simulator read theta directly instead of rebuilding CPD tables.
% --------------------- Pedro Saa UC 2026 ---------------------------------
indexes = getParamIndex(priors);
nNodes  = numel(priors);

off.start = zeros(1,nNodes);
off.rows  = zeros(1,nNodes);
off.cols  = zeros(1,nNodes);

blk = 1;                                                                    % running block counter
for k = 1:nNodes
    off.rows(k)  = size(priors{k},1);                                       % number of parent-state rows
    off.cols(k)  = size(priors{k},2);                                       % number of states of the node
    off.start(k) = indexes(blk,1);
    blk = blk + off.rows(k);
end

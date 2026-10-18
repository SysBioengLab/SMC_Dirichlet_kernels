function S = getHmmSummary(aln,idx)
% Summary statistics of a corpus of pair-HMM alignments (q = 17).
%
% Inputs:   aln  struct from simulateHmmAlignments
%           idx  optional alignment indices, used to bootstrap the corpus
% Output:   S    1 x 17 summary, laid out as
%                  S1 ( 6) emission-symbol counts, symbols 1..6
%                  S2 ( 7) alignment-length histogram: 0, 1, 2, 3-4, 5-8,
%                          9-16, >=17
%                  S3 ( 4) insertion-run-length histogram: 1, 2, 3, >=4
%
% S1 carries only the insertion/match ratio (the emission rows are fixed, so the
% symbol -> state map is deterministic and the chain is symmetric in x- and
% y-insertions). S2 identifies tau on its own. Together they leave one direction
% unidentified: (delta,nu) = (0.02,0.14) and (0.18,0.89) give the same insertion
% fraction and expected length. S3 breaks that ridge, since an insertion run
% continues with probability epsilon = 1-tau-nu, which sweeps 0.76 -> 0.01 along
% it. S3 is deliberately coarser than the transition-count matrix, which would be
% sufficient and would reduce the ABC posterior to the conjugate one.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<2) || isempty(idx)
    L  = aln.L;
    sc = aln.symCounts;
    rc = aln.runCounts;
else
    L  = aln.L(idx);
    sc = aln.symCounts(idx,:);
    rc = aln.runCounts(idx,:);
end

lenEdges = [-0.5 0.5 1.5 2.5 4.5 8.5 16.5 Inf];                            % 0 | 1 | 2 | 3-4 | 5-8 | 9-16 | >=17

S = [sum(sc,1), histcounts(L,lenEdges), sum(rc,1)];

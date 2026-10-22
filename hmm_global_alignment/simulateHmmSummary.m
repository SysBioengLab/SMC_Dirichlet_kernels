function S = simulateHmmSummary(trans,model)
% Simulate a corpus of model.M alignments and summarise it.
%
% Inputs:   trans  4 x 4 transition matrix over [xIns, match, yIns, end]
%           model  struct with fields M, Lmax, p0 and EMIS
% Output:   S      1 x 17 summary (see getHmmSummary)
%
% Wrapper so the observed corpus and every proposed particle follow the same path.
% --------------------- Pedro Saa UC 2026 ---------------------------------
S = getHmmSummary(simulateHmmAlignments(trans,model));

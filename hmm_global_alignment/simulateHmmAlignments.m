function aln = simulateHmmAlignments(trans,model)
% Ancestral sampling of a corpus of pair-HMM alignments.
%
% Inputs:   trans  4 x 4 transition matrix over [xIns, match, yIns, end]
%           model  struct with fields M, Lmax, p0 (1 x 4 initial distribution
%                  over the same states) and EMIS (3 x nSym emission rows of the
%                  three EMITTING states, each summing to one)
% Output:   aln    struct of PER-ALIGNMENT quantities
%                    L          M x 1     number of emitted columns
%                    symCounts  M x nSym  symbol counts
%                    runCounts  M x 4     insertion-run lengths, binned 1,2,3,>=4
%
% Each alignment starts from p0 and emits one column per visit to an emitting
% state until it is absorbed by the end state, so alignments have VARIABLE
% length: L = 0 with probability p0(end), otherwise geometric with rate tau. That
% length distribution is what identifies tau, and it is why the corpus is M
% complete alignments rather than one long sequence.
%
% Written as a direct sampler rather than through hmmgenerate, which has no notion
% of a silent state or of absorption, always produces a fixed number of steps, and
% forces an emission from every state it visits.
%
% All alignments advance together, one column per iteration, so the loop runs as
% many steps as the longest alignment. Lmax caps that length; the truncation is
% applied identically to the observed corpus and to every simulation.
% --------------------- Pedro Saa UC 2026 ---------------------------------
M       = model.M;
Lmax    = model.Lmax;
nSym    = size(model.EMIS,2);
nRunBin = 4;

cumP = cumsum(model.p0);
cumT = cumsum(trans,2);
cumE = cumsum(model.EMIS,2);

s         = 1 + sum(rand(M,1) >= cumP,2);                                  % initial state, 4 = end
L         = zeros(M,1);
symCounts = zeros(M,nSym);
runCounts = zeros(M,nRunBin);
curRun    = zeros(M,1);                                                    % length of the insertion run in progress
alive     = (s<4);

step = 0;
while any(alive) && (step<Lmax)
    idx = find(alive);
    si  = s(idx);

    % Emit one column from the current state
    sym    = 1 + sum(rand(numel(idx),1) >= cumE(si,:),2);
    li     = sub2ind(size(symCounts),idx,sym);                             % one entry per alignment: no repeated indices
    symCounts(li) = symCounts(li) + 1;
    L(idx) = L(idx) + 1;

    % Track insertion runs: a run is a maximal block of consecutive insertion
    % columns, closed by a match column or by the end of the alignment
    ins = (si==1)|(si==3);
    curRun(idx(ins)) = curRun(idx(ins)) + 1;
    cl = idx(~ins);
    cl = cl(curRun(cl)>0);
    if ~isempty(cl)
        lb = sub2ind(size(runCounts),cl,min(curRun(cl),nRunBin));
        runCounts(lb) = runCounts(lb) + 1;
        curRun(cl)    = 0;
    end

    % Transition
    snew   = 1 + sum(rand(numel(idx),1) >= cumT(si,:),2);
    s(idx) = snew;
    fin    = idx(snew==4);
    fin    = fin(curRun(fin)>0);
    if ~isempty(fin)
        lb = sub2ind(size(runCounts),fin,min(curRun(fin),nRunBin));
        runCounts(lb) = runCounts(lb) + 1;
        curRun(fin)   = 0;
    end

    alive(idx) = (snew<4);
    step = step + 1;
end

% Close any run left open by the Lmax truncation
openRuns = find(curRun>0);
if ~isempty(openRuns)
    lb = sub2ind(size(runCounts),openRuns,min(curRun(openRuns),nRunBin));
    runCounts(lb) = runCounts(lb) + 1;
end

aln.L         = L;
aln.symCounts = symCounts;
aln.runCounts = runCounts;

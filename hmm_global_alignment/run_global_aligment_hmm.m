function run_global_aligment_hmm(kernelIdx)
% Run one kernel on the pair-HMM global-alignment benchmark (Fig5.png).
%
% Input:    kernelIdx  1, 3, 4, 6 or 7 (see getKernelSpec)
% Writes:   ./results/results_hmm_full<kernelIdx>.mat, and a checkpoint after
%           every SMC iteration while the run is in progress
%--------------------- Pedro Saa UC 2026 ----------------------------------
clc, clearvars -except kernelIdx

%% ---------------------------- WHAT TO RUN -------------------------------
% Pass the kernel INDEX. One kernel per run; the result file name follows the
% index used in the paper (indices 2 and 5 are unused).
%
%   1  cw_DK          DK          6  beta_logitMGK  BLT-MGK
%   3  MGK            LT-MGK      7  gamma_MGK      LGT-MGK
%   4  ilrt_MGK       ILRT-MGK
%% ------------------------------------------------------------------------

spec = getKernelSpec(kernelIdx);

% Output folder (git does not track empty directories)
if ~isfolder('results'); mkdir('results'); end

% Pool for the SMC proposal loop. Results are independent of pool type and
% worker count (each proposal chunk uses its own RandStream substream). The
% request is clamped to the local profile: 8 workers on the benchmark machine.
setupHmmPool(12);

% Initialize model parameters
rng('default');                                     % for reproducibility
p0     = [.2,.5,.2,.1];                             % initial distribution over [xIns, match, yIns, end]
gam    = .5;                                        % gamma of the paper (gamma() is a builtin)
delta  = .2;
nu     = .8;
epsil  = .1;
tau    = .1;
EMIS   = [.25,zeros(1,5);...                                                % emission probability matrix
          0,.5,.3,.15,.05,0;...
          zeros(1,5),.25];
TRANS  = [epsil,nu,0,tau;delta,gam,delta,tau;...                          % transition probability matrix
          0,nu,epsil,tau;zeros(1,4)];

% The six symbols are alignment-column CATEGORIES, not nucleotides, so the rows
% are normalised (the 0.25 entries encode nucleotide identity, not part of this
% model). test_hmm_predictions scores with the same matrix.
EMIS = EMIS./sum(EMIS,2);

% Training corpus: M complete alignments of VARIABLE length (0 with probability
% p0(end), otherwise geometric with rate tau), which is what identifies tau.
model.M     = 1e3;                                                          % number of alignments
model.Lmax  = 500;                                                          % safety cap on the columns of one alignment
model.p0    = p0;
model.EMIS  = EMIS;

alnData      = simulateHmmAlignments(TRANS,model);
model.xdata  = getHmmSummary(alnData);

% Standard errors of the observed summaries, which weight the discrepancy. No
% closed form applies, so they are bootstrapped by resampling ALIGNMENTS.
rngState = rng;                                                             % keep the bootstrap off the sampler's random stream
rng(20210101);
nBoot = 1000;
Sboot = zeros(nBoot,numel(model.xdata));
for bx = 1:nBoot
    Sboot(bx,:) = getHmmSummary(alnData,randi(model.M,model.M,1));
end
rng(rngState);
model.seData = max(std(Sboot,0,1),1e-8*max(1,abs(model.xdata)));            % guard against a degenerate statistic

% Sampling parameters. rho is in units of standard errors of the q = 17
% summaries, so sqrt(q) ~ 4.1 places every statistic about one standard error
% from its observed value. The floor is simulator noise (~4.9); tolFinal doubles it.
priors{1} = ones(1,2);                                                      % prior definitions
priors{2} = ones(1,3);
alpha     = 50;
N         = 1e4;
pmin_acc  = 1e-5;                                                           % stop if an iteration is provably below this acceptance rate
tolStart  = Inf;
tolFinal  = 12;

%% Run
% The sampler checkpoints every iteration and stops gracefully if acceptance
% collapses, so an incomplete run still leaves a usable population on disk.
checkpointFile = ['./results/checkpoint_',spec.resultFile,'.mat'];
resultPath     = ['./results/',spec.resultFile,'.mat'];

fprintf('\nKernel %d: %s (%s), M = %g alignments, N = %g, tolFinal = %g\n',...
        spec.idx,spec.kernel,spec.paperName,model.M,N,tolFinal);

try
    population = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,...
                                             spec.kernel,spec.method,checkpointFile,pmin_acc);
catch err
    fprintf(2,'\nRun failed: %s\n%s\n',err.identifier,err.message);
    if isfile(checkpointFile)
        fprintf(2,'The last completed iteration is preserved in %s\n',checkpointFile);
    end
    rethrow(err);
end

% Save under the variable name expected by the post-processing scripts
out.(spec.varName) = population;
save(resultPath,'-struct','out');

fprintf('Saved %s (converged = %d, final tolerance = %g, ESS = %.1f, time = %.1f s)\n',...
        resultPath,population.converged,population.finalTol,population.ESS(end),population.time);
end

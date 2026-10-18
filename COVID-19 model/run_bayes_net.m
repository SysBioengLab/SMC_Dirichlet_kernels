function run_bayes_net(kernelIdx)
% Run one kernel on the COVID-19 screening Bayes net benchmark (Fig6.png).
%
% Input:    kernelIdx  1, 3, 4, 6 or 7 (see getKernelSpec)
% Writes:   ./results/results_bayes_net<kernelIdx>.mat, and a checkpoint after
%           every SMC iteration while the run is in progress
%------------------------- Pedro Saa UC 2026 ------------------------------
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
setupCovidPool(12);

% Initialize parameters
rng('default');                                                              % for reproducibility
load('xdata.mat','xdata');                                                   % n x 14 cohort: COV, 7 symptoms, 6 odour outcomes

% Read the data
xdata = xdata+1;                                                             % convention false==1, true==2
nSize = size(xdata,1);

% Observed summaries and their standard errors, which weight the discrepancy.
% Bootstrapped rather than binomial, since the last four statistics are central
% moments; for the proportions the bootstrap reproduces the binomial value.
Sdata = getCovidSummary(xdata);

rngState = rng;                                                              % keep the bootstrap off the sampler's random stream
rng(20210101);
nBoot = 1000;
Sboot = zeros(nBoot,numel(Sdata));
for bx = 1:nBoot
    Sboot(bx,:) = getCovidSummary(xdata(randi(nSize,nSize,1),:));
end
rng(rngState);

seData = std(Sboot,0,1);
seData = max(seData,1e-8*max(1,abs(Sdata)));                                 % guard against a degenerate statistic

% Initialize sampling parameters
priors_sim = initializePriors_final_sim(xdata);

% Sampling parameters. rho is in units of standard errors of the q = 31
% summaries. Two floors bound it: simulator noise (~7.1) and model-data mismatch
% (~22), so anything below ~23 is unreachable; prior draws give rho ~80. N must
% support a 37- to 72-dimensional kernel covariance (gamma_MGK needs N_alive > 72).
alpha     = 50;
N         = 1e4;
pmin_acc  = 1e-5;                                                            % stop if an iteration is provably below this acceptance rate
tolStart  = Inf;
tolFinal  = 35;

%% Run
% The sampler checkpoints every iteration and stops gracefully if acceptance
% collapses, so an incomplete run still leaves a usable population on disk.
checkpointFile = ['./results/checkpoint_',spec.resultFile,'.mat'];
resultPath     = ['./results/',spec.resultFile,'.mat'];

fprintf('\nKernel %d: %s (%s), N = %g, tolFinal = %g\n',...
        spec.idx,spec.kernel,spec.paperName,N,tolFinal);

try
    population = bayes_net_example(priors_sim,N,alpha,tolStart,tolFinal,...
                                   spec.kernel,spec.method,Sdata,seData,nSize,checkpointFile,pmin_acc);
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

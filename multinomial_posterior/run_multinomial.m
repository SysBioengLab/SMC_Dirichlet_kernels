function run_multinomial(kernelIdx)
% Run one kernel on the multinomial benchmark with a known posterior (Fig4.png).
%
% Input:    kernelIdx  1, 3, 4, 6 or 7 (see getKernelSpec)
% Writes:   ./results/results_multinomial<kernelIdx>.mat, and a checkpoint after
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

% This example runs SERIALLY: a proposal costs microseconds, so parfor dispatch
% overhead dominates. Reproducibility comes from the seed below.
rng('default');

% Six-component multinomial. The counts are sufficient, so the posterior is
% conjugate: theta ~ Dir(prior), y | theta ~ Multi(n,theta) => theta | y ~
% Dir(prior + y). Here Dir(ones(1,6)) and y = [1 2 3 4 5 5] give Dir([2 3 4 5 6 6]).
priorAlpha  = ones(1,6);
model.xdata = [1,2,3,4,5,5];
model.n     = sum(model.xdata);
model.posteriorAlpha = priorAlpha + model.xdata;                            % the reference every kernel is measured against

% Initialize sampling parameters
priors{1} = priorAlpha;                                                     % a single Dirichlet block
alpha     = 50;
N         = 1e4;
pmin_acc  = 1e-5;                                                           % stop if an iteration is provably below this acceptance rate
tolStart  = Inf;

% rho is the l1 norm of the count differences, an even integer in [0,2n]. With
% tolFinal = 0 only exact reproductions of the observed counts are accepted, so
% the ABC posterior is the exact posterior rather than an approximation of it.
tolFinal  = 0;

%% Run
% The sampler checkpoints every iteration and stops gracefully if acceptance
% collapses, so an incomplete run still leaves a usable population on disk.
checkpointFile = ['./results/checkpoint_',spec.resultFile,'.mat'];
resultPath     = ['./results/',spec.resultFile,'.mat'];

fprintf('\nKernel %d: %s (%s), n = %g trials, N = %g, tolFinal = %g\n',...
        spec.idx,spec.kernel,spec.paperName,model.n,N,tolFinal);
fprintf('Exact posterior: Dir(%s)\n',num2str(model.posteriorAlpha));

try
    population = multinomial_example(model,priors,N,alpha,tolStart,tolFinal,...
                                     spec.kernel,spec.method,checkpointFile,pmin_acc);
catch err
    fprintf(2,'\nRun failed: %s\n%s\n',err.identifier,err.message);
    if isfile(checkpointFile)
        fprintf(2,'The last completed iteration is preserved in %s\n',checkpointFile);
    end
    rethrow(err);
end

% Carry the reference posterior with the results
population.posteriorAlpha = model.posteriorAlpha;
population.priorAlpha     = priorAlpha;

out.(spec.varName) = population;
save(resultPath,'-struct','out');

fprintf('Saved %s (converged = %d, final tolerance = %g, ESS = %.1f, time = %.1f s)\n',...
        resultPath,population.converged,population.finalTol,population.ESS(end),population.time);
end

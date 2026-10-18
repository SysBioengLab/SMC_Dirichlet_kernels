function M = getPopulationMetrics(population,label)
% Performance metrics of one SMC-ABC run, as reported in the paper.
%
% Input:    population  struct returned by any of the ASMC samplers
%           label       optional display name (defaults to the kernel name)
% Output:   M           scalar struct, one field per metric
%
%   AAR    average acceptance rate over the SMC iterations. Entries equal to 1
%          are excluded: they belong to the prior-sampling stage, where every
%          draw is accepted when tolStart = Inf.
%   ESS    effective sample size of the final population, 1/sum(w^2)
%   time   seconds spent in the SMC stage
%   perf   time/ESS, the paper's cputime/ESS shown in the computational-
%          performance panel of every benchmark figure.
%          Lower is better; it penalises a kernel that buys acceptance by
%          proposing locally (high AAR, low particle diversity).
%
% Particles live in SMCparams (Bayes net, pair-HMM) or thetaSMC (multinomial).
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<2) || isempty(label)
    if isfield(population,'kernel')
        label = population.kernel;
    else
        label = sprintf('method %d',population.method);
    end
end

if isfield(population,'SMCparams')
    theta = population.SMCparams;
else
    theta = population.thetaSMC;
end

w   = population.weights{end};
w   = w/sum(w);
ESS = 1/sum(w.^2);

aar = population.acceptRate(:)';
aar(aar==1) = [];                                                          % drop the prior stage (and any iteration with nothing to replenish)
if isempty(aar); aar = NaN; end

M.label       = label;
M.nParticles  = size(theta,1);
M.nParams     = size(theta,2);
M.nPopulations = numel(population.tolSMC);                                 % prior population + SMC iterations
M.nSMCiter    = numel(population.tolSMC)-1;
M.AAR         = mean(aar);
M.minAcceptRate = min(aar);
M.ESS         = ESS;
M.ESSfraction = ESS/M.nParticles;
M.time        = population.time;
M.perf        = population.time/ESS;                                       % cputime/ESS
M.tolStart    = population.tolSMC(1);
M.tolFinal    = population.tolSMC(end);

if isfield(population,'converged')
    M.converged = population.converged;
else
    M.converged = NaN;                                                     % sampler predates the flag
end
if isfield(population,'shrinkDelta')
    d = population.shrinkDelta(:)';
    M.shrinkFraction = mean(d>0);                                          % share of iterations needing covariance shrinkage
    M.shrinkMax      = max(d);
else
    M.shrinkFraction = NaN;
    M.shrinkMax      = NaN;
end

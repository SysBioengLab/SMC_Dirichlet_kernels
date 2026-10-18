function population = multinomial_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method,checkpointFile,pmin_acc)
% SMC-ABC for the 6-component multinomial model with a generic kernel choice.
%
% Inputs:   model     struct with fields n (number of trials) and xdata
%                     (observed category counts, which sum to n)
%           priors    cell of Dirichlet hyper-parameters, one per block
%           N         number of particles
%           alpha     percentage of "alive" particles retained per iteration
%           tolStart  initial tolerance (use Inf)
%           tolFinal  target tolerance
%           kernel    'cw_DK' | 'MGK' | 'ilrt_MGK' | 'beta_logitMGK' |
%                     'gamma_MGK'
%           method    DK tuning method (1 MLE, 2 min-KL, 3 total variance);
%                     pass [] for the MGK kernels
%           checkpointFile  optional .mat path; when given, the population is
%                     written there after every completed SMC iteration
%           pmin_acc  optional minimum acceptance rate (default 1e-5). The run
%                     stops gracefully once an iteration is provably below it
%
% Output:   population  scalar struct; its fields are listed in scripts/README.md
%
% Each proposal is scored by drawing model.n multinomial trials and comparing the
% resulting counts with the observed ones through the l1 norm. The counts are
% sufficient for theta, so with tolFinal = 0 the ABC posterior is the exact
% posterior Dir(prior + xdata) and any departure from it is a defect of the
% sampler rather than of the approximation.
%--------------------- Pedro Saa UC 2026 ----------------------------------
% Validated up front: an unrecognised name would fall through to the Gaussian
% branch and fail later with a misleading "not positive semi-definite".
validKernels = {'cw_DK','MGK','ilrt_MGK','beta_logitMGK','gamma_MGK'};
if ~any(strcmp(kernel,validKernels))
    error('multinomial_example:unknownKernel',...
          'Unknown kernel ''%s''. Valid kernels are: %s.',kernel,strjoin(validKernels,', '));
end

if (nargin<9); checkpointFile = ''; end
if (nargin<10) || isempty(pmin_acc); pmin_acc = 1e-5; end

% Initializing sampling parameters
currTol   = tolStart;                                                      % initial tolerance
numPriors = numel(priors);

% Flat-vector layout of the Dirichlet blocks, and of their images in the
% transformed space (one coordinate fewer per block)
indexes   = zeros(numPriors,2);
transIdxs = zeros(numPriors,2);
posFull   = 1;
posTrans  = 1;
for ix = 1:numPriors
    dBlock          = numel(priors{ix});
    indexes(ix,:)   = [posFull,posFull+dBlock-1];
    transIdxs(ix,:) = [posTrans,posTrans+dBlock-2];
    posFull         = posFull + dBlock;
    posTrans        = posTrans + dBlock - 1;
end
D = indexes(end);

if (numel(model.xdata)~=indexes(end)) || (sum(model.xdata)~=model.n)
    error('multinomial_example:summarySize',...
         ['model.xdata must hold one count per category (%d expected) and sum to '...
          'model.n = %d; it has %d entries summing to %d.'],...
          indexes(end),model.n,numel(model.xdata),sum(model.xdata));
end

rho       = zeros(N,1);
weights   = ones(N,1)/N;
SMCparams = zeros(N,D);

% Placeholders so the per-iteration state can be packed uniformly for every
% kernel; each kernel fills in only the fields it actually uses.
mappingFxn = []; invMapping = []; invlogitMapping = []; logitMapping = [];
gammaMapping = []; log_pz = []; newPriors = {}; fullMapping = [];
alphaDK = {}; logBetaDK = {}; cholK = []; logNormK = [];

% Set up the space each MGK variant works in: mappings to and from the
% simplex, and the prior density expressed in the transformed coordinates
if strcmp(kernel,'gamma_MGK')

    % Initialize array with log( SMCparams )
    SMCparamsLog = zeros(N,D);

    % Construct prior distributions
    for ix = 1:numPriors
        newPriors{ix} = priors{ix};

        % Write the mapping and inverse mapping functions
        if (ix==1)
            invMapping   = @(y) exp(y)/sum(exp(y));
            gammaMapping = @(x) log(x);

            % Log-Gamma density of z = log(v), v ~ Gamma(k,1):
            %   log p(z) = k*z - exp(z) - log(Gamma(k))
            log_pz = @(z,k) k(:).*z(:) - exp(z(:)) - gammaln(k(:));
        end
    end

elseif strcmp(kernel,'beta_logitMGK')

    % Construct beta (D-1)-dimensional priors
    for ix = 1:numPriors
        numParams     = numel(priors{ix});
        newPriors{ix} = zeros(numParams-1,2);
        for jx = 1:numParams-1
            newPriors{ix}(jx,:) = [sum(priors{ix}(jx+1:end)),priors{ix}(jx)];
        end

        % Build inverse mapping function (from ind. betas to Dirichlet on the simplex)
        if (numParams>2)
            invMapping{ix} = @(y) [1-y(:,1),(exp(cumsum(log(y(:,1:end-1)),2))'.*(1-y(:,2:end))')',prod(y,2)];
        elseif (numParams==2)
            invMapping{ix} = @(y) [1-y(:,1),prod(y,2)];
        else
            invMapping{ix} = @(y) (1-y(:,1));
        end

        % Logit mapping and the resulting log-prior density, built once
        if (ix==1)
            invlogitMapping = @(z) exp(z)./(1+exp(z));
            logitMapping    = @(y) log(y./(1-y));

            % Logit-Beta density of z = logit(y), y ~ Beta(a,b):
            %   log p(z) = a*z - (a+b)*log(1+exp(z)) - log B(a,b)
            log_pz = @(z,alpha) alpha(:,1).*z(:) - sum(alpha,2).*log(1 + exp(z(:))) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));
        end
    end

elseif strcmp(kernel,'MGK')
    U    = [];
    Uinv = [];
    b    = [];
    for ix = 1:numPriors
        numParams = numel(priors{ix});
        U    = blkdiag(U,[eye(numParams-1),zeros(numParams-1,1)]);
        Uinv = blkdiag(Uinv,[eye(numParams-1);-ones(1,numParams-1)]);
        b    = [b;zeros(numParams-1,1);1];
    end
    mappingFxn = @(x) (U*x')';
    invMapping = @(y) (Uinv*y' + b(:,ones(size(y,1),1)))';

elseif strcmp(kernel,'ilrt_MGK')

    % Build mapping function to the unconstrained space
    U_full{numPriors}     = [];
    mappingFxn{numPriors} = [];
    invMapping{numPriors} = [];
    detJ{numPriors}       = [];
    log_pz{numPriors}     = [];
    for ix = 1:numPriors
        numParams = numel(priors{ix});
        U_ilrt    = [];
        for jx = 1:numParams-1
            U_ilrt(:,jx) = sqrt(jx/(jx+1))*[(1/jx)*ones(jx,1);-1;zeros(numParams-jx-1,1)];
        end
        U_full{ix}  = [U_ilrt,ones(numParams,1)/sqrt(numParams)];
        fullMapping = blkdiag(fullMapping,U_ilrt);

        % Build mapping functions
        mappingFxn{ix} = @(x) log(x)*U_ilrt;
        invMapping{ix} = @(y) (exp(U_ilrt*y'))'./repmat(sum((exp(U_ilrt*y'))',2),1,numParams);

        % Formulate determinant of the transformation's Jacobian for the new density (we map z onto the full space)
        detJ{ix} = @(z) det(U_full{ix}.*exp(U_full{ix}*repmat((log(invMapping{ix}(z))*U_full{ix})',1,numParams)));

        % Transformed log-prior density (z mapped onto the full space).
        % NOTE: detJ is that of the FULL d x d map, while the d-1 free coordinates
        % need the (d-1)-dimensional one; they differ by log(sqrt(d)) per block.
        % Constant in z, so it cancels on normalisation: log_pz is proportional
        % to the transformed prior.
        log_pz{ix} = @(z,alpha) (alpha-1)*(U_full{ix}*(log(invMapping{ix}(z))*U_full{ix})') + log(detJ{ix}(z)) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));
    end
end

%% I: Prior sampling
popIdx     = 1;
counter    = 1;
totalCount = 0;
while (counter<=N)

    % Update total count
    totalCount = totalCount + 1;

    % Check type of kernel for sampling from the prior
    if strcmp(kernel,'ilrt_MGK')

        % Exact draw from the prior: in z-space it is the pushforward of
        % Dir(alpha) under z = log(x)*U_ilrt, so sampling the Dirichlet and
        % mapping it is i.i.d. -- no Markov chain needed.
        z_curr{numPriors} = [];
        for ix = 1:numPriors
            x          = randg(priors{ix});
            x          = x/sum(x);
            z_curr{ix} = mappingFxn{ix}(x);
        end

        % Map the draw onto the simplex
        theta = zeros(1,D);
        for ix = 1:numPriors
            theta(indexes(ix,1):indexes(ix,2)) = invMapping{ix}(z_curr{ix});
        end

    elseif strcmp(kernel,'beta_logitMGK')||strcmp(kernel,'gamma_MGK')

        % Exact draws from the prior: both are products of standard univariate
        % distributions composed with a bijection, so no Markov chain is needed.
        %   beta_logitMGK : y_j ~ Beta(alphabar_j,alpha_j), then z = logit(y)
        %   gamma_MGK     : v_j ~ Gamma(alpha_j,1),         then z = log(v)
        z_curr{numPriors} = [];
        for ix = 1:numPriors
            if strcmp(kernel,'beta_logitMGK')
                g          = randg(newPriors{ix});                         % newPriors{ix} is (D-1)x2 = [alphabar,alpha]
                ybeta      = g(:,1)./sum(g,2);                             % Beta variate per component
                z_curr{ix} = logitMapping(ybeta)';
            else
                z_curr{ix} = gammaMapping(randg(newPriors{ix}));           % log of independent Gamma variates
            end
        end

        % Map the draw onto the simplex
        theta = zeros(1,D);
        if strcmp(kernel,'beta_logitMGK')
            for ix = 1:numPriors
                ybeta = invlogitMapping(z_curr{ix});
                theta(indexes(ix,1):indexes(ix,2)) = invMapping{ix}(ybeta);
            end
        else
            thetaLog = zeros(1,D);
            for ix = 1:numPriors
                thetaLog(indexes(ix,1):indexes(ix,2)) = z_curr{ix};
                theta(indexes(ix,1):indexes(ix,2))    = invMapping(z_curr{ix});
            end
        end

    else
        % All other kernels: draw from the original Dirichlet prior
        theta = zeros(1,D);
        for ix = 1:numPriors
            probTheta = randg(priors{ix});
            theta(indexes(ix,1):indexes(ix,2)) = probTheta/sum(probTheta);
        end
    end

    % Simulate a multinomial dataset from the model and summarise it
    simulated = simulateMultinomialSummary(theta,model);

    % Accept/Reject particle (l1 norm of the count differences)
    rhoNew = sum(abs(simulated-model.xdata));
    if (rhoNew<=currTol)
        rho(counter)         = rhoNew;
        SMCparams(counter,:) = theta;
        if strcmp(kernel,'gamma_MGK')
            SMCparamsLog(counter,:) = thetaLog;
        end
        population.rho{popIdx}(counter) = rhoNew;
        counter = counter + 1;
    end
end

% Sort particles by discrepancy
[rho,xorder] = sort(rho,'ascend');

% Store the prior population
population.weights{1}     = weights(xorder);
population.acceptRate(1)  = (counter-1)/totalCount;
population.ESS(1)         = N;
population.rho{1}         = rho;
population.Kparams{1}     = nan;
population.shrinkDelta(1) = 0;                                             % no kernel is tuned during the prior phase
population.rejectBudget(1,:) = [totalCount,0,totalCount-(counter-1),counter-1];   % prior draws are always inside the simplex
population.tolSMC(1)      = max(rho);
SMCparams                 = SMCparams(xorder,:);
if strcmp(kernel,'gamma_MGK')
    SMCparamsLog          = SMCparamsLog(xorder,:);
end

%% II: SMC step
% Main loop
tstartSMC = tic;
while (currTol~=tolFinal)&&(population.acceptRate(end)>pmin_acc)

    % Advance one iteration
    popIdx = popIdx + 1;

    % The percentile normally drives the schedule; the additive term binds only
    % if it fails to move, guaranteeing a decreasing sequence and termination.
    % The step is 2, not 1, because rho lives on a lattice of spacing 2: it is the
    % l1 norm of the difference between two count vectors that both sum to
    % model.n, so it is always EVEN. Odd tolerances are unattainable, and an
    % iteration targeting one finds nothing to replenish.
    currTol = max([min([currTol-2,prctile(population.rho{popIdx-1},alpha)]),tolFinal]);

    % Snapshot of the last COMPLETE particle set. If this iteration cannot be
    % filled, the partial replenishment is discarded and this set is returned.
    % Only the particles need saving: the stored weights, discrepancies and
    % schedule of the previous iteration are already in population.
    prevSMCparams = SMCparams;

    % Drop particles below the tolerance
    population.tolSMC(popIdx)     = currTol;
    droppedParticles              = find(rho>currTol);
    rho(droppedParticles)         = [];
    weights(droppedParticles)     = [];
    SMCparams(droppedParticles,:) = [];
    weights                       = weights/sum(weights);
    if strcmp(kernel,'gamma_MGK')
        SMCparamsLog(droppedParticles,:) = [];
    end

    % Start of next sampling step
    Nreplenish       = N-numel(weights);                                                            % Particles to replenish
    Nalive           = N-Nreplenish;
    weightsReplenish = zeros(Nreplenish,1);
    SMCreplenish     = zeros(Nreplenish,D);
    rhoReplenish     = zeros(Nreplenish,1);
    if strcmp(kernel,'gamma_MGK')                                                                   % only this kernel carries a log-space copy
        SMCreplenishLog = zeros(Nreplenish,D);
    end

    % Compute kernel parameters for the joint
    if strcmp(kernel,'gamma_MGK')                                                                   % work on the log-parameters
        SMCparams_original = SMCparams;
        SMCparams = SMCparamsLog;

    elseif strcmp(kernel,'MGK')                                                                     % work only with free parameters
        SMCparams_original = SMCparams;
        SMCparams = mappingFxn(SMCparams);

    elseif strcmp(kernel,'ilrt_MGK')                                                                % work only with free parameters on the transformed space
        SMCparams_original = SMCparams;
        SMCparams = log(SMCparams)*fullMapping;

    elseif strcmp(kernel,'beta_logitMGK')
        SMCparams_original = SMCparams;
        SMCtemp = [];
        for ix = 1:numPriors
            SMCtemp = [SMCtemp,mappingFxn_beta(SMCparams(:,indexes(ix,1):indexes(ix,2)))];
        end
        SMCparams = logitMapping(SMCtemp);
        clearvars('SMCtemp');
    end

    % Compute kernel parameters
    [Kparam,shrinkDelta] = getGeneralKernelParam(SMCparams,priors,kernel,method,weights,alpha);
    population.Kparams{popIdx}     = Kparam;
    population.shrinkDelta(popIdx) = shrinkDelta;                                                   % 0 unless block-diagonal shrinkage was needed

    % Weight terms depending only on the alive population and the kernel
    % hyper-parameter are constant over the replenishment loop, so they are
    % hoisted here rather than rebuilt for every accepted particle.
    logW = log(weights);
    if strcmp(kernel,'cw_DK')
        alphaDK   = cell(numPriors,1);
        logBetaDK = cell(numPriors,1);
        for ix = 1:numPriors
            alphaDK{ix}   = Kparam(ix)*SMCparams(:,indexes(ix,1):indexes(ix,2))+1;
            logBetaDK{ix} = sum(gammaln(alphaDK{ix}),2) - gammaln(sum(alphaDK{ix},2));
        end
    else
        % Gaussian kernels: factorise the covariance once and work in logs.
        % cholcov(.,0) is what mvnpdf itself uses; unlike chol it tolerates the
        % rank-deficient covariance a concentrated population produces.
        [cholK,errFlag] = cholcov(Kparam,0);
        if (errFlag~=0)
            error('multinomial_example:kernelCovariance',...
                  'Kernel covariance is not positive semi-definite (kernel %s, iteration %d).',kernel,popIdx);
        end
        logNormK = -0.5*size(Kparam,1)*log(2*pi) - sum(log(diag(cholK)));
    end

    % Sampling loop, run SERIALLY: a single proposal costs microseconds, so parfor
    % dispatch overhead dominates (measured ~47% core utilisation and a longer
    % wall time). Reproducibility comes from the global stream seeded in the driver.
    cumw = cumsum([0,weights']);                                                                % hoisted out of the proposal loop

    P             = struct();
    P.kernel      = kernel;
    P.indexes     = indexes;
    P.transIdxs   = transIdxs;
    P.numPriors   = numPriors;
    P.priors      = priors;
    P.newPriors   = newPriors;
    P.Kparam      = Kparam;
    P.logW        = logW;
    P.SMCparams   = SMCparams;
    P.cumw        = cumw;
    P.model       = model;
    P.currTol     = currTol;
    P.alphaDK     = alphaDK;
    P.logBetaDK   = logBetaDK;
    P.cholK       = cholK;
    P.logNormK    = logNormK;
    P.mappingFxn  = mappingFxn;
    P.invMapping  = invMapping;
    P.invlogit    = invlogitMapping;
    P.logitMapping = logitMapping;
    P.log_pz      = log_pz;

    nAcc       = 0;
    totalCount = 0;
    budget     = zeros(1,3);                                                                    % [infeasible, tolRejected, accepted] over this iteration
    collapsed  = false;

    while (nAcc<Nreplenish)
        [ok,th,thLog,rh,lw,outcome] = proposeParticle(P);
        totalCount    = totalCount + 1;
        budget(outcome) = budget(outcome) + 1;

        if ok
            nAcc = nAcc + 1;
            SMCreplenish(nAcc,:)   = th;
            rhoReplenish(nAcc)     = rh;
            weightsReplenish(nAcc) = lw;
            if strcmp(kernel,'gamma_MGK')
                SMCreplenishLog(nAcc,:) = thLog;
            end
        end

        % Stop gracefully once the iteration is provably below pmin_acc, judged
        % only after 1e6 proposals so the rate is meaningful. There is no cap on
        % total effort: a slow run above pmin_acc is allowed to finish. The
        % partial batch is discarded rather than mixed in (see the rollback).
        if (totalCount>=1e6) && (nAcc/totalCount<pmin_acc)
            warning('multinomial_example:acceptanceCollapse',...
                    ['Acceptance collapsed (kernel %s, iteration %d): %d proposals, %d accepted '...
                     '(rate %.2e < pmin %.0e). Returning the population from iteration %d.'],...
                     kernel,popIdx,totalCount,nAcc,nAcc/totalCount,pmin_acc,popIdx-1);
            collapsed = true;
            break
        end
    end

    if collapsed
        % Roll back to the last complete population and leave the main loop.
        % population.weights/rho for this iteration are written only after this
        % point, so they still describe the previous one; restore the particles
        % to match them and drop this iteration from the schedule.
        SMCparams = prevSMCparams;
        population.tolSMC(popIdx)      = [];
        population.Kparams(popIdx)     = [];
        population.shrinkDelta(popIdx) = [];
        popIdx  = popIdx-1;
        currTol = population.tolSMC(popIdx);
        break
    end

    % Re-normalize importance weights for robust calculations
    weightsReplenish = weightsReplenish - max(weightsReplenish);
    weightsReplenish = exp(weightsReplenish);

    % Normalize the weights over all the particles
    weightsReplenish = weightsReplenish/sum(weightsReplenish);
    weights          = [Nalive*weights;Nreplenish*weightsReplenish]/N;
    rho              = [rho;rhoReplenish];

    % Append new parameters according to their size
    if strcmp(kernel,'cw_DK')
        SMCparams = [SMCparams;SMCreplenish];
    elseif strcmp(kernel,'gamma_MGK')
        SMCparamsLog = [SMCparamsLog;SMCreplenishLog];
        SMCparams    = [SMCparams_original;SMCreplenish];
    else
        SMCparams = [SMCparams_original;SMCreplenish];
    end

    % Sort particles by discrepancy
    [rho,xorder] = sort(rho,'ascend');
    SMCparams    = SMCparams(xorder,:);
    weights      = weights(xorder);
    if strcmp(kernel,'gamma_MGK')
        SMCparamsLog = SMCparamsLog(xorder,:);
    end

    % Store this iteration
    population.weights{popIdx}    = weights;
    population.acceptRate(popIdx) = Nreplenish/totalCount;
    population.ESS(popIdx)        = 1/sum(weights.^2);
    population.rho{popIdx}        = rho;
    population.rejectBudget(popIdx,:) = [totalCount,budget];                                    % [proposals, infeasible, tolRejected, accepted]
    if isnan(population.acceptRate(popIdx)); population.acceptRate(popIdx) = 1; end             % No particle needed replenishing (Nreplenish = 0), so nothing was rejected

    % Checkpoint after every completed iteration: a run that is interrupted (or
    % killed) still leaves a usable population on disk.
    if ~isempty(checkpointFile)
        writeCheckpoint(checkpointFile,population,SMCparams,...
                        kernel,method,currTol,tolFinal,toc(tstartSMC),numel(model.xdata));
    end
end

% Collect final samples
population.time   = toc(tstartSMC);
population.method = method;
population.kernel = kernel;

% The population is returned whether or not the target tolerance was reached, so
% check converged first: an early stop samples a looser ABC posterior.
population.converged = (currTol==tolFinal);
population.finalTol  = currTol;
population.SMCparams = SMCparams;
population.nStats    = numel(model.xdata);                                 % stamps the summary design, so runs from different setups are not mixed

if ~population.converged
    disp(['Acceptance rate below minimum, kernel: ',kernel,'. ',...
          'Target tolerance NOT reached (final tolerance = ',num2str(currTol),...
          ', target = ',num2str(tolFinal),'). Population returned with converged = false.']);
end

end

%% ------------------------------------------------------------------------
function writeCheckpoint(fileName,population,SMCparams,kernel,method,currTol,tolFinal,elapsed,nStats)
% Write the population obtained so far, in the same shape the driver saves at
% the end, so a checkpoint can be analysed exactly like a completed run.
population.SMCparams = SMCparams;
population.kernel    = kernel;
population.method    = method;
population.time      = elapsed;
population.finalTol  = currTol;
population.converged = (currTol==tolFinal);
population.nStats    = nStats;
save(fileName,'population');
end

%% ------------------------------------------------------------------------
function lp = logMvn(X,mu,cholK,logNormK)
% log N(X_i ; mu, K) per row of X, from a Cholesky factor computed once per SMC
% iteration. Equals log(mvnpdf(X,mu,K)) without refactorising K on every call.
Z  = (X - mu)/cholK;
lp = logNormK - 0.5*sum(Z.^2,2);
end

%% ------------------------------------------------------------------------
function [ok,th,thLog,rh,lw,outcome] = proposeParticle(P)
% Generate ONE proposal, returning: acceptance flag, the particle, its log-space
% copy (gamma_MGK only), its discrepancy, its un-normalised log importance weight
% and the outcome code.
%
% outcome is 1, 2 or 3 for the three mutually exclusive fates of a proposal --
% outside the simplex, rejected by the tolerance, accepted -- which separates a
% kernel that proposes inadmissible moves from one that merely proposes far.
nP    = P.numPriors;
D     = P.indexes(end);
ok    = false;
th    = zeros(1,D);
thLog = zeros(1,D);
rh    = 0;
lw    = 0;
newProbTheta = cell(nP,1);

% Resample from the previous population
prevIdx      = sum(rand(1)>=P.cumw);
prevParticle = P.SMCparams(prevIdx,:);

% Perform perturbation according to the kernel
thetaProposal    = zeros(1,D);
thetaProposalLog = zeros(1,D);
newParticle      = [];
if strcmp(P.kernel,'cw_DK')                                                                 % DK perturbation on each Dir param
    for ix = 1:nP
        probProposal = randg(P.Kparam(ix)*prevParticle(P.indexes(ix,1):P.indexes(ix,2))+1);
        probProposal = probProposal/sum(probProposal);
        thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = probProposal;
        newProbTheta{ix} = probProposal;
    end

elseif strcmp(P.kernel,'MGK')                                                               % MGK on the Dir composition
    newParticle = mvnrnd(prevParticle,P.Kparam);
    newParticle = P.invMapping(newParticle);
    for ix = 1:nP
        newProbTheta{ix} = newParticle(P.indexes(ix,1):P.indexes(ix,2));
        thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = newProbTheta{ix};
    end

elseif strcmp(P.kernel,'ilrt_MGK')||strcmp(P.kernel,'beta_logitMGK')
    newParticle = mvnrnd(prevParticle,P.Kparam);                                            % MGK on the transformed space
    for ix = 1:nP
        if ~strcmp(P.kernel,'beta_logitMGK')
            newProbTheta{ix} = P.invMapping{ix}(newParticle(P.transIdxs(ix,1):P.transIdxs(ix,2)));
        else
            ybeta = P.invlogit(newParticle(P.transIdxs(ix,1):P.transIdxs(ix,2)));
            newProbTheta{ix} = P.invMapping{ix}(ybeta);
        end
        thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = newProbTheta{ix};
    end

elseif strcmp(P.kernel,'gamma_MGK')
    newParticle      = mvnrnd(prevParticle,P.Kparam);
    thetaProposalLog = newParticle;
    for ix = 1:nP
        newProbTheta{ix} = P.invMapping(newParticle(P.indexes(ix,1):P.indexes(ix,2)));
        thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = newProbTheta{ix};
    end
end

% Checkpoint for feasibility of the proposal
if any(thetaProposal<0); outcome = 1; return; end

% Simulate a multinomial dataset from the model and summarise it
simulated = simulateMultinomialSummary(thetaProposal,P.model);
rhoNew    = sum(abs(simulated-P.model.xdata));
if (rhoNew>P.currTol); outcome = 2; return; end

% Compute the (log)weight for the accepted particle
if strcmp(P.kernel,'cw_DK')||strcmp(P.kernel,'MGK')
    probPrior = 0;
    for ix = 1:nP
        probPrior = probPrior + probDirichlet(P.priors{ix},newProbTheta{ix},1);
    end
    if strcmp(P.kernel,'cw_DK')
        particleWeight = P.logW;
        for ix = 1:nP
            particleWeight = particleWeight + (P.alphaDK{ix}-1)*log(newProbTheta{ix}(:)) - P.logBetaDK{ix};
        end
        particleWeight = logsumexp(particleWeight);
    else
        particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,P.mappingFxn(newParticle),P.cholK,P.logNormK));
    end

elseif strcmp(P.kernel,'ilrt_MGK')
    probPrior = 0;
    for ix = 1:nP
        probPrior = probPrior + P.log_pz{ix}(P.mappingFxn{ix}(newProbTheta{ix}),P.priors{ix});
    end
    particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,newParticle,P.cholK,P.logNormK));

elseif strcmp(P.kernel,'beta_logitMGK')
    % newParticle IS the z-space draw the kernel produced, so the prior is
    % evaluated there directly rather than through an identity round trip.
    probPrior = 0;
    for ix = 1:nP
        zBlock    = newParticle(P.transIdxs(ix,1):P.transIdxs(ix,2));
        probPrior = probPrior + sum(P.log_pz(zBlock,P.newPriors{ix}));
    end
    particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,newParticle,P.cholK,P.logNormK));

else                                                                                        % gamma_MGK
    probPrior = 0;
    for ix = 1:nP
        zBlock    = thetaProposalLog(P.indexes(ix,1):P.indexes(ix,2));
        probPrior = probPrior + sum(P.log_pz(zBlock,P.newPriors{ix}));
    end
    particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,thetaProposalLog,P.cholK,P.logNormK));
end

outcome = 3;
ok      = true;
th      = thetaProposal;
thLog   = thetaProposalLog;
rh      = rhoNew;
lw      = probPrior - particleWeight;
end

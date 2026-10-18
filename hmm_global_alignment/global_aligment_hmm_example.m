function population = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method,checkpointFile,pmin_acc)
% SMC-ABC for the pair-HMM global alignment model with a generic kernel choice.
%
% Inputs:   model     struct with fields M (number of alignments), Lmax, p0,
%                     EMIS (normalised emission rows of the emitting states),
%                     xdata (observed summary) and seData (its standard errors)
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
% Each proposed particle is scored by simulating a corpus of model.M complete
% alignments, summarising it (q = 17: symbol counts, alignment-length histogram
% and insertion-run-length histogram) and comparing with the observed summary
% through a standard-error weighted Euclidean distance, so the discrepancy is in
% units of standard errors.
%--------------------- Pedro Saa UC 2026 ----------------------------------
% Validated up front: an unrecognised name would fall through to the Gaussian
% branch and fail later with a misleading "not positive semi-definite".
validKernels = {'cw_DK','MGK','ilrt_MGK','beta_logitMGK','gamma_MGK'};
if ~any(strcmp(kernel,validKernels))
    error('global_aligment_hmm_example:unknownKernel',...
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

if ~isfield(model,'seData') || (numel(model.seData)~=numel(model.xdata))
    error('global_aligment_hmm_example:summarySize',...
         ['model.xdata (%d entries) and model.seData must be the observed summary and '...
          'its standard errors, both produced by getHmmSummary. Regenerate them in the '...
          'driver.'],numel(model.xdata));
end

rho       = zeros(N,1);
weights   = ones(N,1)/N;
SMCparams = zeros(N,D);

% Define auxiliary functions
% x(1) = gamma , x(2) = delta , x(3) = nu , x(4) = epsil , x(5) = tau
TRANS    = @(x) [x(4),x(3),0,x(5);x(2),x(1),x(2),x(5);0,x(3),x(4),x(5);zeros(1,4)];
hmm_pars = @(theta1,theta2) [theta1(:,1).*(1-theta2(:,3)),theta1(:,2).*(1-theta2(:,3))/2,theta2];

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

        % Transformed log-prior density (z mapped onto the full space). detJ is
        % that of the FULL d x d map; the d-1 free coordinates need the
        % (d-1)-dimensional one. They differ by log(sqrt(d)) per block, constant
        % in z, so log_pz is proportional to the transformed prior.
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

        % Exact i.i.d. draw: in z-space the prior is the pushforward of
        % Dir(alpha) under z = log(x)*U_ilrt, so a Dirichlet draw suffices.
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

        % Exact draws: both priors are products of standard univariate
        % distributions composed with a bijection.
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

    % Simulate a corpus of alignments from the model and summarise it
    hmm_theta = hmm_pars(theta(indexes(1,1):indexes(1,2)),theta(indexes(2,1):indexes(2,2)));
    simulated = simulateHmmSummary(TRANS(hmm_theta),model);

    % Accept/Reject particle (standard-error weighted Euclidean discrepancy)
    rhoNew = sqrt(sum(((simulated-model.xdata)./model.seData).^2));
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
tolStep   = 1e-3;                                                          % minimum RELATIVE decrease of the tolerance per iteration
tstartSMC = tic;
while (currTol~=tolFinal)&&(population.acceptRate(end)>pmin_acc)

    % Advance one iteration
    popIdx = popIdx + 1;

    % The percentile normally drives the schedule; the multiplicative term binds
    % only if it fails to move, guaranteeing a decreasing sequence and
    % termination. Multiplicative rather than additive because rho is now in
    % units of standard errors (single digits), not raw symbol counts.
    currTol = max([min([currTol*(1-tolStep),prctile(population.rho{popIdx-1},alpha)]),tolFinal]);

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

    % Terms depending only on the alive population and the kernel
    % hyper-parameter are constant over the replenishment loop, so hoist them.
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
            error('global_aligment_hmm_example:kernelCovariance',...
                  'Kernel covariance is not positive semi-definite (kernel %s, iteration %d).',kernel,popIdx);
        end
        logNormK = -0.5*size(Kparam,1)*log(2*pi) - sum(log(diag(cholK)));
    end

    % Sampling loop. Proposals within an iteration are mutually independent, so
    % the loop is embarrassingly parallel. Work is issued in CHUNKS to build one
    % RandStream per parfor task rather than per proposal; each chunk uses its own
    % deterministic substream, making a run reproducible and INDEPENDENT OF THE
    % WORKER COUNT. With no pool open, parfor runs serially.
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
    P.TRANS       = TRANS;
    P.hmm_pars    = hmm_pars;
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

    chunkSize  = 32;
    pHat       = max(population.acceptRate(popIdx-1),1e-3);                                     % acceptance estimate from the previous iteration
    nAcc       = 0;
    totalCount = 0;
    budget     = zeros(1,3);                                                                    % [infeasible, tolRejected, accepted] over this iteration
    collapsed  = false;
    subBase    = popIdx*1e6;                                                                    % distinct substream block per SMC iteration

    while (nAcc<Nreplenish)
        need   = Nreplenish-nAcc;
        nChunk = max(1,ceil(1.5*need/(pHat*chunkSize)));
        nChunk = min(nChunk,4096);                                                              % cap the batch so memory stays bounded

        Ok = cell(nChunk,1); Th = cell(nChunk,1); TL = cell(nChunk,1);
        Rh = cell(nChunk,1); Lw = cell(nChunk,1); Cn = cell(nChunk,1);
        parfor c = 1:nChunk
            [Ok{c},Th{c},TL{c},Rh{c},Lw{c},Cn{c}] = proposeChunk(P,chunkSize,subBase+c);
        end
        budget = budget + sum(vertcat(Cn{:}),1);
        okB = vertcat(Ok{:}); thB = vertcat(Th{:}); tlB = vertcat(TL{:});
        rhB = vertcat(Rh{:}); lwB = vertcat(Lw{:});

        totalCount = totalCount + nChunk*chunkSize;
        subBase    = subBase + nChunk;

        sel  = find(okB);
        take = min(numel(sel),Nreplenish-nAcc);
        if (take>0)
            keepIdx = sel(1:take);
            SMCreplenish(nAcc+1:nAcc+take,:)   = thB(keepIdx,:);
            rhoReplenish(nAcc+1:nAcc+take)     = rhB(keepIdx);
            weightsReplenish(nAcc+1:nAcc+take) = lwB(keepIdx);
            if strcmp(kernel,'gamma_MGK')
                SMCreplenishLog(nAcc+1:nAcc+take,:) = tlB(keepIdx,:);
            end
            nAcc = nAcc + take;
        end
        pHat = max(numel(sel)/(nChunk*chunkSize),1e-4);                                         % refresh the estimate for the next batch

        % Stop gracefully once the iteration is provably below pmin_acc, judged
        % only after 1e6 proposals so the rate is meaningful. There is no cap on
        % total effort: a slow run above pmin_acc is allowed to finish. The
        % partial batch is discarded rather than mixed in (see the rollback).
        if (totalCount>=1e6) && (nAcc/totalCount<pmin_acc)
            warning('global_aligment_hmm_example:acceptanceCollapse',...
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

    % Checkpoint every iteration, so an interrupted run still leaves a
    % usable population on disk.
    if ~isempty(checkpointFile)
        writeCheckpoint(checkpointFile,population,SMCparams,hmm_pars,indexes,...
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
population.nStats    = numel(model.xdata);                                 % stamps the summary/discrepancy setup, so post-processing can refuse to mix runs from different designs

% HMM transition parameters are a deterministic function of the particles, so
% they are derived once here rather than tracked through the sampler.
population.HMMparams = hmm_pars(SMCparams(:,indexes(1,1):indexes(1,2)),SMCparams(:,indexes(2,1):indexes(2,2)));

if ~population.converged
    disp(['Acceptance rate below minimum, kernel: ',kernel,'. ',...
          'Target tolerance NOT reached (final tolerance = ',num2str(currTol),...
          ', target = ',num2str(tolFinal),'). Population returned with converged = false.']);
end

end

%% ------------------------------------------------------------------------
function writeCheckpoint(fileName,population,SMCparams,hmm_pars,indexes,kernel,method,currTol,tolFinal,elapsed,nStats)
% Write the population obtained so far in the same shape the driver saves at the
% end, so a checkpoint can be analysed exactly like a completed run.
population.SMCparams = SMCparams;
population.HMMparams = hmm_pars(SMCparams(:,indexes(1,1):indexes(1,2)),SMCparams(:,indexes(2,1):indexes(2,2)));
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
function [ok,th,thLog,rh,lw,cnt] = proposeChunk(P,m,subIdx)
% Generate m independent proposals from one deterministic substream, returning
% per proposal: acceptance flag, particle, log-space copy (gamma_MGK only),
% discrepancy and un-normalised log importance weight. All inputs are fixed for
% the SMC iteration, so chunks may run in any order; subIdx selects the substream
% that makes the result independent of how they are scheduled.
%
% cnt splits the proposals into the three mutually exclusive outcomes [outside the
% simplex, rejected by the tolerance, accepted], which separates a kernel that
% proposes inadmissible moves from one that merely proposes far.
st = RandStream('mrg32k3a','Seed',20210101);
st.Substream = subIdx;
RandStream.setGlobalStream(st);

nP = P.numPriors;
D  = P.indexes(end);
ok = false(m,1); th = zeros(m,D); thLog = zeros(m,D); rh = zeros(m,1); lw = zeros(m,1);
cnt = zeros(1,3);                                                                                % [infeasible, tolRejected, accepted]
newProbTheta = cell(nP,1);

for q = 1:m

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
    if any(thetaProposal<0); cnt(1) = cnt(1)+1; continue; end

    % Simulate a corpus of alignments from the model and summarise it
    hmm_theta = P.hmm_pars(thetaProposal(P.indexes(1,1):P.indexes(1,2)),thetaProposal(P.indexes(2,1):P.indexes(2,2)));
    simulated = simulateHmmSummary(P.TRANS(hmm_theta),P.model);
    rhoNew    = sqrt(sum(((simulated-P.model.xdata)./P.model.seData).^2));
    if (rhoNew>P.currTol); cnt(2) = cnt(2)+1; continue; end

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

    cnt(3)     = cnt(3)+1;
    ok(q)      = true;
    th(q,:)    = thetaProposal;
    thLog(q,:) = thetaProposalLog;
    rh(q)      = rhoNew;
    lw(q)      = probPrior - particleWeight;
end
end

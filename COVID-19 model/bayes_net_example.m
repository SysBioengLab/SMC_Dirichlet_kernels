function population = bayes_net_example(priors,N,alpha,tolStart,tolFinal,kernel,method,Sdata,seData,nSubj,checkpointFile,pmin_acc)
% SMC-ABC for the COVID-19 screening Bayes net with a generic kernel choice.
%
% Inputs:   priors    15x1 cell of Dirichlet hyper-parameters (counts)
%           N         number of particles
%           alpha     percentage of "alive" particles retained per iteration
%           tolStart  initial tolerance (use Inf)
%           tolFinal  target tolerance, in units of standard errors
%           kernel    'cw_DK' | 'MGK' | 'ilrt_MGK' | 'beta_logitMGK' |
%                     'gamma_MGK'
%           method    DK tuning method (1 MLE, 2 min-KL, 3 total variance);
%                     pass [] for the MGK kernels
%           Sdata     1x31 summary statistics of the observed cohort
%           seData    1x31 standard errors of the observed statistics
%           nSubj     cohort size to simulate (use the size of the data)
%           checkpointFile  optional .mat path; when given, the population is
%                     written there after every completed SMC iteration
%           pmin_acc  optional minimum acceptance rate (default 1e-5). The run
%                     stops gracefully once an iteration is provably below it
%
% Output:   population  scalar struct; its fields are listed in scripts/README.md
%
% Each proposal is scored by simulating a synthetic cohort of nSubj subjects
% (ancestral sampling), summarising it and comparing with the observed summaries
% through a standard-error weighted Euclidean distance. Discrepancies are thus in
% units of standard errors, so sqrt(31) ~ 5.6 places every statistic about one
% standard error from its observed value.
% --------------------- Pedro Saa UC 2026 ---------------------------------
% Validated up front: an unrecognised name would fall through to the Gaussian
% branch and fail later with a misleading "not positive semi-definite".
validKernels = {'cw_DK','MGK','ilrt_MGK','beta_logitMGK','gamma_MGK'};
if ~any(strcmp(kernel,validKernels))
    error('bayes_net_example:unknownKernel',...
          'Unknown kernel ''%s''. Valid kernels are: %s.',kernel,strjoin(validKernels,', '));
end

if (nargin<11); checkpointFile = ''; end
if (nargin<12) || isempty(pmin_acc); pmin_acc = 1e-5; end

% Initializing sampling parameters
currTol         = tolStart;                                                % initial tolerance
rho             = zeros(N,1);
weights         = ones(N,1)/N;
indexes         = getParamIndex(priors);
transIdxs       = getTransformedIndex(indexes);
numDirParams    = size(indexes,1);
SMCparams       = zeros(N,indexes(end));
samplingPriors  = reshapeProbs(priors,priors,numDirParams,'sampling');
covidOffsets    = getCovidOffsets(priors);                                 % flat-vector CPD offsets, used by simulateCovidSummary

% Placeholders so the per-iteration state can be packed uniformly for every
% kernel; each kernel fills in only the fields it actually uses.
mappingFxn = []; invMapping = []; invlogitMapping = []; logitMapping = [];
gammaMapping = []; log_pz = []; newPriors = {};
alphaDK = {}; logBetaDK = {}; cholK = []; logNormK = [];

% Set up the space each MGK variant works in: mappings to and from the
% simplex, and the prior density expressed in the transformed coordinates
if strcmp(kernel,'gamma_MGK')

    % Log-space copy of the particles, the coordinates gamma_MGK perturbs
    SMCparamsLog = zeros(N,indexes(end));

    % Construct prior distributions
    for ix = 1:numDirParams
        newPriors{ix} = samplingPriors{ix};

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
    for ix = 1:numDirParams
        numParams     = numel(samplingPriors{ix});
        newPriors{ix} = zeros(numParams-1,2);
        for jx = 1:numParams-1
            newPriors{ix}(jx,:) = [sum(samplingPriors{ix}(jx+1:end)),samplingPriors{ix}(jx)];
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
    for ix = 1:numDirParams
        numParams = numel(samplingPriors{ix});
        U    = blkdiag(U,[eye(numParams-1),zeros(numParams-1,1)]);
        Uinv = blkdiag(Uinv,[eye(numParams-1);-ones(1,numParams-1)]);
        b    = [b;zeros(numParams-1,1);1];
    end
    mappingFxn = @(x) (U*x')';
    invMapping = @(y) (Uinv*y' + b(:,ones(size(y,1),1)))';

elseif strcmp(kernel,'ilrt_MGK')

    % Build mapping function to the unconstrained space
    fullMapping              = [];
    U_full{numDirParams}     = [];
    mappingFxn{numDirParams} = [];
    invMapping{numDirParams} = [];
    z_curr{numDirParams}     = [];
    detJ{numDirParams}       = [];
    log_pz{numDirParams}     = [];
    for ix = 1:numDirParams
        numParams = numel(samplingPriors{ix});
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
        z_curr{numDirParams} = [];
        for ix = 1:numDirParams
            x          = randg(samplingPriors{ix});
            x          = x/sum(x);
            z_curr{ix} = mappingFxn{ix}(x);
        end

        % Map the draw onto the simplex
        theta = zeros(1,indexes(end));
        probTheta{numDirParams} = [];
        for ix = 1:numDirParams
            probTheta{ix} = invMapping{ix}(z_curr{ix});
            theta(indexes(ix,1):indexes(ix,2)) = probTheta{ix};
        end

    elseif strcmp(kernel,'beta_logitMGK')||strcmp(kernel,'gamma_MGK')

        % Exact draws: both priors are products of standard univariate
        % distributions composed with a bijection.
        %   beta_logitMGK : y_j ~ Beta(alphabar_j,alpha_j), then z = logit(y)
        %   gamma_MGK     : v_j ~ Gamma(alpha_j,1),         then z = log(v)
        z_curr{numDirParams} = [];
        for ix = 1:numDirParams
            if strcmp(kernel,'beta_logitMGK')
                g          = randg(newPriors{ix});                                  % newPriors{ix} is (D-1)x2 = [alphabar,alpha]
                ybeta      = g(:,1)./sum(g,2);                                      % Beta variate per component
                z_curr{ix} = logitMapping(ybeta)';
            elseif strcmp(kernel,'gamma_MGK')
                z_curr{ix} = gammaMapping(randg(newPriors{ix}));                    % log of independent Gamma variates
            end
        end

        % Map the draw onto the simplex
        theta = zeros(1,indexes(end));
        probTheta{numDirParams} = [];
        if strcmp(kernel,'beta_logitMGK')
            for ix = 1:numDirParams
                ybeta         = invlogitMapping(z_curr{ix});
                probTheta{ix} = invMapping{ix}(ybeta);
                theta(indexes(ix,1):indexes(ix,2)) = probTheta{ix};
            end

        elseif strcmp(kernel,'gamma_MGK')
            thetaLog = zeros(1,indexes(end));
            for ix = 1:numDirParams
                thetaLog(indexes(ix,1):indexes(ix,2)) = z_curr{ix};
                probTheta{ix} = invMapping(z_curr{ix});
                theta(indexes(ix,1):indexes(ix,2)) = probTheta{ix};
            end
        end

    else
        % All other kernels: draw from the original Dirichlet prior
        theta = zeros(1,indexes(end));
        probTheta{numDirParams} = [];
        for ix = 1:numDirParams
            probTheta{ix} = randg(samplingPriors{ix});
            probTheta{ix} = probTheta{ix}/sum(probTheta{ix});
            theta(indexes(ix,1):indexes(ix,2)) = probTheta{ix};
        end
    end

    % Simulate a synthetic cohort from the model and summarise it in one pass,
    % reading the flat parameter vector directly (no CPD cell array is built).
    Ssim = simulateCovidSummary(theta,covidOffsets,nSubj);

    % Accept/Reject particle (standard-error weighted Euclidean discrepancy)
    rhoNew = sqrt(sum(((Ssim-Sdata)./seData).^2));
    if (rhoNew<=currTol)
        rho(counter)          = rhoNew;
        SMCparams(counter,:)  = theta;
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
population.weights{1}    = weights(xorder);
population.acceptRate(1) = (counter-1)/totalCount;
population.ESS(1)        = N;
population.rho{1}        = rho;
population.Kparams{1}    = nan;
population.shrinkDelta(1) = 0;                                                                      % no kernel is tuned during the prior phase
population.rejectBudget(1,:) = [totalCount,0,totalCount-(counter-1),counter-1];                     % prior draws are always inside the simplex
population.tolSMC(1)     = max(rho);
SMCparams                = SMCparams(xorder,:);
if strcmp(kernel,'gamma_MGK')
    SMCparamsLog         = SMCparamsLog(xorder,:);
end

%% II: SMC step
% Main loop
tolStep   = 1e-3;                                                                                  % minimum RELATIVE decrease of the tolerance per iteration
tstartSMC = tic;
while (currTol~=tolFinal)&&(population.acceptRate(end)>pmin_acc)

    % Advance one iteration
    popIdx = popIdx + 1;

    % The percentile normally drives the schedule; the multiplicative term binds
    % only if it fails to move, guaranteeing a decreasing sequence and termination.
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
    SMCreplenish     = zeros(Nreplenish,indexes(end));
    rhoReplenish     = zeros(Nreplenish,1);
    if strcmp(kernel,'gamma_MGK')                                                                   % only this kernel carries a log-space copy
        SMCreplenishLog = zeros(Nreplenish,indexes(end));
    end

    % Compute kernel parameters for the joint
    if strcmp(kernel,'gamma_MGK')                                                                   % work on the log-transformed parameters
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
        for ix = 1:numDirParams
            SMCtemp = [SMCtemp,mappingFxn_beta(SMCparams(:,indexes(ix,1):indexes(ix,2)))];
        end
        SMCparams = logitMapping(SMCtemp);
        clearvars('SMCtemp');
    end

    % Compute kernel parameters
    [Kparam,shrinkDelta] = getGeneralKernelParam(SMCparams,kernel,method,weights,alpha,indexes,numDirParams);
    population.Kparams{popIdx}    = Kparam;
    population.shrinkDelta(popIdx) = shrinkDelta;                                                   % 0 unless block-diagonal shrinkage was needed

    % Terms depending only on the alive population and the kernel
    % hyper-parameter are constant over the replenishment loop, so hoist them.
    logW = log(weights);
    if strcmp(kernel,'cw_DK')
        alphaDK   = cell(numDirParams,1);
        logBetaDK = cell(numDirParams,1);
        for ix = 1:numDirParams
            alphaDK{ix}   = Kparam(ix)*SMCparams(:,indexes(ix,1):indexes(ix,2))+1;
            logBetaDK{ix} = sum(gammaln(alphaDK{ix}),2) - gammaln(sum(alphaDK{ix},2));
        end
    else
        % Gaussian kernels: factorise the covariance once and work in logs, to
        % avoid mvnpdf underflowing to 0 in 37-72 dimensions. cholcov(.,0) is what
        % mvnpdf itself uses; unlike chol it tolerates the rank-deficient
        % covariance a concentrated population produces.
        [cholK,errFlag] = cholcov(Kparam,0);
        if (errFlag~=0)
            error('bayes_net_example:kernelCovariance',...
                  'Kernel covariance is not positive semi-definite (kernel %s, iteration %d).',kernel,popIdx);
        end
        logNormK = -0.5*size(Kparam,1)*log(2*pi) - sum(log(diag(cholK)));
    end

    % Sampling loop. Proposals within an iteration are mutually independent, so
    % the loop is embarrassingly parallel. Work is issued in CHUNKS to build one
    % RandStream per parfor task rather than per proposal (~40 us against ~300 us
    % of useful work); each chunk uses its own deterministic substream, making a
    % run reproducible and INDEPENDENT OF THE WORKER COUNT. With no pool open,
    % parfor degrades to a serial loop.
    cumw = cumsum([0,weights']);                                                                % hoisted out of the proposal loop

    P                = struct();
    P.kernel         = kernel;
    P.indexes        = indexes;
    P.transIdxs      = transIdxs;
    P.numDirParams   = numDirParams;
    P.samplingPriors = samplingPriors;
    P.newPriors      = newPriors;
    P.Kparam         = Kparam;
    P.logW           = logW;
    P.SMCparams      = SMCparams;
    P.cumw           = cumw;
    P.covidOffsets   = covidOffsets;
    P.nSubj          = nSubj;
    P.Sdata          = Sdata;
    P.seData         = seData;
    P.currTol        = currTol;
    P.alphaDK        = alphaDK;
    P.logBetaDK      = logBetaDK;
    P.cholK          = cholK;
    P.logNormK       = logNormK;
    P.mappingFxn     = mappingFxn;
    P.invMapping     = invMapping;
    P.invlogit       = invlogitMapping;
    P.logitMapping   = logitMapping;
    P.log_pz         = log_pz;

    chunkSize  = 64;
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
            warning('bayes_net_example:acceptanceCollapse',...
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
    if strcmp(kernel,'MGK')||strcmp(kernel,'ilrt_MGK')||strcmp(kernel,'beta_logitMGK')
        SMCparams = [SMCparams_original;SMCreplenish];
    elseif strcmp(kernel,'gamma_MGK')
        SMCparamsLog = [SMCparamsLog;SMCreplenishLog];
        SMCparams    = [SMCparams_original;SMCreplenish];
    else
        SMCparams = [SMCparams;SMCreplenish];
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
    population.rejectBudget(popIdx,:) = [totalCount,budget];                                      % [proposals, infeasible, tolRejected, accepted]
    if isnan(population.acceptRate(popIdx)); population.acceptRate(popIdx) = 1; end                % No particle needed replenishing (Nreplenish = 0), so nothing was rejected

    % Checkpoint every iteration, so an interrupted run still leaves a
    % usable population on disk.
    if ~isempty(checkpointFile)
        writeCheckpoint(checkpointFile,population,SMCparams,kernel,method,currTol,tolFinal,toc(tstartSMC));
    end
end

% Collect final samples
population.time   = toc(tstartSMC);
population.method = method;
population.kernel = kernel;

% The population is returned whether or not the target tolerance was reached, so
% check converged first: an early stop samples a looser ABC posterior.
population.converged = (currTol==tolFinal);
population.finalTol   = currTol;
population.SMCparams = SMCparams;

if ~population.converged
    disp(['Acceptance rate below minimum, kernel: ',kernel,'. ',...
          'Target tolerance NOT reached (final tolerance = ',num2str(currTol),...
          ', target = ',num2str(tolFinal),'). Population returned with converged = false.']);
end

end

%% ------------------------------------------------------------------------
function writeCheckpoint(fileName,population,SMCparams,kernel,method,currTol,tolFinal,elapsed)
% Write the population obtained so far in the same shape the driver saves at the
% end, so a checkpoint can be analysed exactly like a completed run.
population.SMCparams = SMCparams;
population.kernel    = kernel;
population.method    = method;
population.time      = elapsed;
population.finalTol  = currTol;
population.converged = (currTol==tolFinal);
save(fileName,'population');
end

%% ------------------------------------------------------------------------
function lp = logMvn(X,mu,cholK,logNormK)
% log N(X_i ; mu, K) per row of X, from a Cholesky factor computed once per SMC
% iteration. Equals log(mvnpdf(X,mu,K)) without refactorising K on every call or
% underflowing to -Inf in high dimensions.
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

nD = P.numDirParams;
D  = P.indexes(end);
ok = false(m,1); th = zeros(m,D); thLog = zeros(m,D); rh = zeros(m,1); lw = zeros(m,1);
cnt = zeros(1,3);                                                                                % [infeasible, tolRejected, accepted]
newProbTheta = cell(nD,1);

for q = 1:m

    % Resample from the previous population
    prevIdx      = sum(rand(1)>=P.cumw);
    prevParticle = P.SMCparams(prevIdx,:);

    % Perform perturbation according to the kernel
    thetaProposal    = zeros(1,D);
    thetaProposalLog = zeros(1,D);
    newParticle      = [];
    if strcmp(P.kernel,'cw_DK')                                                                 % DK perturbation on each Dir param
        for ix = 1:nD
            probProposal = randg(P.Kparam(ix)*prevParticle(P.indexes(ix,1):P.indexes(ix,2))+1);
            probProposal = probProposal/sum(probProposal);
            thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = probProposal;
            newProbTheta{ix} = probProposal;
        end

    elseif strcmp(P.kernel,'MGK')                                                               % MGK on the Dir composition
        newParticle = mvnrnd(prevParticle,P.Kparam);
        newParticle = P.invMapping(newParticle);
        for ix = 1:nD
            newProbTheta{ix} = newParticle(P.indexes(ix,1):P.indexes(ix,2));
            thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = newProbTheta{ix};
        end

    elseif strcmp(P.kernel,'ilrt_MGK')||strcmp(P.kernel,'beta_logitMGK')
        newParticle = mvnrnd(prevParticle,P.Kparam);                                            % MGK on the transformed space
        for ix = 1:nD
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
        for ix = 1:nD
            newProbTheta{ix} = P.invMapping(newParticle(P.indexes(ix,1):P.indexes(ix,2)));
            thetaProposal(P.indexes(ix,1):P.indexes(ix,2)) = newProbTheta{ix};
        end
    end

    % Checkpoint for feasibility of the proposal
    if any(thetaProposal<0); cnt(1) = cnt(1)+1; continue; end

    % Simulate a synthetic cohort from the model and summarise it in one pass
    Ssim   = simulateCovidSummary(thetaProposal,P.covidOffsets,P.nSubj);
    rhoNew = sqrt(sum(((Ssim-P.Sdata)./P.seData).^2));
    if (rhoNew>P.currTol); cnt(2) = cnt(2)+1; continue; end

    % Compute the (log)weight for the accepted particle
    if strcmp(P.kernel,'cw_DK')||strcmp(P.kernel,'MGK')
        probPrior = 0;
        for ix = 1:nD
            probPrior = probPrior + probDirichlet(P.samplingPriors{ix},newProbTheta{ix},1);
        end
        if strcmp(P.kernel,'cw_DK')
            particleWeight = P.logW;
            for ix = 1:nD
                particleWeight = particleWeight + (P.alphaDK{ix}-1)*log(newProbTheta{ix}(:)) - P.logBetaDK{ix};
            end
            particleWeight = logsumexp(particleWeight);
        else
            particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,P.mappingFxn(newParticle),P.cholK,P.logNormK));
        end

    elseif strcmp(P.kernel,'ilrt_MGK')
        probPrior = 0;
        for ix = 1:nD
            probPrior = probPrior + P.log_pz{ix}(P.mappingFxn{ix}(newProbTheta{ix}),P.samplingPriors{ix});
        end
        particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,newParticle,P.cholK,P.logNormK));

    elseif strcmp(P.kernel,'beta_logitMGK')
        % newParticle IS the z-space draw the kernel produced, so the prior is
        % evaluated there directly rather than through an identity round trip.
        probPrior = 0;
        for ix = 1:nD
            zBlock    = newParticle(P.transIdxs(ix,1):P.transIdxs(ix,2));
            probPrior = probPrior + sum(P.log_pz(zBlock,P.newPriors{ix}));
        end
        particleWeight = logsumexp(P.logW + logMvn(P.SMCparams,newParticle,P.cholK,P.logNormK));

    else                                                                                        % gamma_MGK
        probPrior = 0;
        for ix = 1:nD
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

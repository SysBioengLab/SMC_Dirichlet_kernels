function population = bayes_net_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method)
% This scripts runs the single-reaction model
% Inputs:
%
% Outputs: population (structure array)
% --------------------- Pedro Saa UC 2021 ---------------------------------
% Initializing sampling parameters
pmin_acc        = 1e-5;                                                    % min. acceptance probability
currTol         = tolStart;                                                % initial tolerance
rho             = zeros(N,1);
ESS             = N;
weights         = ones(N,1)/N;
indexes         = getParamIndex(priors);
transIdxs       = getTransformedIndex(indexes);
numDirParams    = size(indexes,1);
HMM_params{N,1} = [];
SMCparams       = zeros(N,indexes(end));
samplingPriors  = reshapeProbs(priors,priors,numDirParams,'sampling');

% Initalize special transformation features for particular kernels
if strcmp(kernel,'gamma_MGK')
    
% Initialize array with log( SMCparams )
    SMCparamsLog = zeros(N,indexes(end));
    
    % Construct prior distributions    
    for ix = 1:numDirParams
        newPriors{ix} = samplingPriors{ix};        
        
        % Write mapping and inverse maping functions
        if (ix==1)
            invMapping   = @(y) exp(y)/sum(exp(y));
            gammaMapping = @(x) log(x);
            
            % Formulate transformed log-prior density
            %               k*z  - exp(z) - log( gamma( k ) )
            log_pz = @(z,k) k(:).*z(:) - exp(z(:)) - gammaln(k(:));
            
            % Initialize parameters for the Markov Chain (sampling directly from the prior is not possible in this case)
            burnin   = 5e3;
            thinning = 1e2;
        end
        
        % Generate a suitable covariance matrix by sampling from the
        % original prior, mapping the variables onto the transformed space,
        % and then computing the sample covariance matrix
        x          = randg(newPriors{ix}(ones(1e4,1),:));
        y          = gammaMapping(x);
        Q{ix}      = cov(y);                                                 % transtion kernel MVN(0,sigma) (symmetric)
        z_curr{ix} = zeros(1,numel(newPriors{ix}));%mean(y);        
    end
elseif strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')
    
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
        
                % If we are using the logit transformation, formulate logit
        % tranformation fnxs
        if strcmp(kernel,'beta_logitMGK')
            
            % Write mapping and inverse maping functions
            if (ix==1)
                invlogitMapping = @(z) exp(z)./(1+exp(z));
                logitMapping    = @(y) log(y./(1-y));
                
                % Formulate transformed log-prior density (we map z onto the full space)
                %                        alpha*z     -    (alpha+beta)*log( )           -          log ( betaFxn )
                log_pz = @(z,alpha) alpha(:,1).*z(:) - sum(alpha,2).*log(1 + exp(z(:))) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));
                
                % Initialize parameters for the Markov Chain (sampling directly from the prior is not possible in this case)
                burnin   = 2e3;
                thinning = 5e1;
            end
            
            % Generate a suitable covariance matrix by sampling from the
            % original prior, mapping the variables onto the transformed space,
            % and then computing the sample covariance matrix
            Q{ix}      = [];
            z_curr{ix} = [];
            for jx = 1:size(newPriors{ix},1)
                y = randg(repmat(newPriors{ix}(jx,:),1e4,1));
                y = y(:,1)./sum(y,2);
                z = logitMapping(y);
                Q{ix} = [Q{ix};cov(z)];                                            % transtion kernel N(0,sigma) (symmetric)
                z_curr{ix} = [z_curr{ix},mean(z)];
            end
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
    Q{numDirParams}          = [];
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
        
        % Formulate transformed log-prior density (we map z onto the full space)
        %               log( exp(             sum( )       )         )                      +  log( det_term ) -     log ( betaFxn )
        log_pz{ix} = @(z,alpha) (alpha-1)*(U_full{ix}*(log(invMapping{ix}(z))*U_full{ix})') + log(detJ{ix}(z)) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));
        
        % Initialize state of the Markov Chain
        z_curr{ix} = mappingFxn{ix}(samplingPriors{ix}/sum(samplingPriors{ix}));
        
        % Generate a suitable covariance matrix by sampling from the
        % original prior, mapping the variables onto the transformed space,
        % and then computing the sample covariance matrix
        x = randg(samplingPriors{ix}(ones(1e4,1),:));
        x = x./repmat(sum(x,2),1,numel(samplingPriors{ix}));
        z = mappingFxn{ix}(x);
        Q{ix} = cov(z);                                 % transtion kernel MGK(0,sigma) (symmetric)
    end
    
    % Initialize parameters for the Markov Chain (sampling directly from the prior is not possible in this case)
    burnin   = 2e3;
    thinning = 5e1;
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
        
        % Start MCMC from the previous point and check if burn-in is
        % required
        if (burnin>0)
            numSteps = burnin;
            burnin   = 0;                              % No more burn-in
        else
            numSteps = thinning;
        end
        
        % Run Metropolis-Hastings as necessary        
        for ix = 1:numSteps
            z_cand{numDirParams} = [];
            log_pz_cand = 0;
            log_pz_curr = 0;
            for jx = 1:numDirParams
                z_cand{jx} = mvnrnd(z_curr{jx},Q{jx});
                log_pz_cand = log_pz_cand + log_pz{jx}(z_cand{jx},samplingPriors{jx});
                log_pz_curr = log_pz_curr + log_pz{jx}(z_curr{jx},samplingPriors{jx});
            end
                        
            % Compute acceptance prob
            accProb = exp( log_pz_cand - log_pz_curr);
            
            % Accept proposal
            if (rand(1) <= accProb)
                z_curr = z_cand;                
            end
        end
        
        % Map proposals to the simplex
        theta = [];
        probTheta{numDirParams} = [];
        for ix = 1:numDirParams
            probTheta{ix} = invMapping{ix}(z_curr{ix});
            theta = [theta,probTheta{ix}];
        end
        
        % Sample from beta/logit or log-gamma prior if required
    elseif strcmp(kernel,'beta_logitMGK')||strcmp(kernel,'gamma_MGK')
        
        % Start MCMC from the previous point and check if burn-in is
        % required
        if (burnin>0)
            numSteps = burnin;
            burnin   = 0;                              % No more burn-in
        else
            numSteps = thinning;
        end
        
        % Run Metropolis-Hastings
        for ix = 1:numSteps            
            z_cand{numDirParams} = [];
            log_pz_cand = 0;
            log_pz_curr = 0;
            for jx = 1:numDirParams
                if strcmp(kernel,'beta_logitMGK')
                    z_cand{jx} = mvnrnd(z_curr{jx},diag(Q{jx}));
                elseif strcmp(kernel,'gamma_MGK')
                    z_cand{jx} = mvnrnd(z_curr{jx},Q{jx});
                end
                log_pz_cand = log_pz_cand + sum(log_pz(z_cand{jx},newPriors{jx}));
                log_pz_curr = log_pz_curr + sum(log_pz(z_curr{jx},newPriors{jx}));
            end

            % Compute acceptance prob
            accProb = exp( log_pz_cand - log_pz_curr );
            
            % Accept proposal
            if (rand(1) <= accProb)
                z_curr = z_cand;
            end
        end
        
        % Map proposals to the simplex
        theta = [];
        probTheta{numDirParams} = [];
        if strcmp(kernel,'beta_logitMGK')
            for ix = 1:numDirParams
                ybeta         = invlogitMapping(z_curr{ix});
                probTheta{ix} = invMapping{ix}(ybeta);
                theta         = [theta,probTheta{ix}];
            end

        elseif strcmp(kernel,'gamma_MGK')
            thetaLog = [];
            for ix = 1:numDirParams
                thetaLog      = [thetaLog,z_curr{ix}];
                probTheta{ix} = invMapping(z_curr{ix});
                theta         = [theta,probTheta{ix}];
            end
        end

        % Sample from beta prior if required
    elseif strcmp(kernel,'beta_MGK')
        
        % Draw from (D-1) independent beta distributions for each prior
        theta = [];
        probTheta{numDirParams} = [];
        for ix = 1:numDirParams
            temp = randg(newPriors{ix});           
            temp = (temp(:,1)./sum(temp,2))';
            probTheta{ix} = invMapping{ix}(temp);
            theta         = [theta,probTheta{ix}];
        end

        % Sample from the original prior for remaining cases
    else
        theta = [];
        probTheta{numDirParams} = [];
        for ix = 1:numDirParams
            probTheta{ix} = randg(samplingPriors{ix});
            probTheta{ix} = probTheta{ix}/sum(probTheta{ix});
            theta         = [theta,probTheta{ix}];
        end
    end

    % Calculate HMM parameters
    simProbTheta = reshapeProbs(probTheta,priors,numDirParams,'simulation');
    simulated    = getSummaryStatistics(simulateBayesNet(simProbTheta,model.n));

    % Accept/Reject particle
    rhoNew = sqrt(sum((simulated(:)-model.Sdata(:)).^2)/sum(model.Sdata(:)~=0));
    if (rhoNew<=currTol)
        rho(counter)          = rhoNew;
        SMCparams(counter,:)  = theta;
        if strcmp(kernel,'gamma_MGK')
            SMCparamsLog(counter,:) = thetaLog;
        end
		HMM_params{counter,1} = probTheta;
        population.rho{popIdx}(counter) = rhoNew;
        counter = counter + 1;
    end
end

% Sort particles by tolerance score
[rho,xorder] = sort(rho,'ascend');

% Make final assignation
population.weights{1}    = weights(xorder);
population.acceptRate(1) = (counter-1)/totalCount;
population.ESS(1)        = ESS;
population.rho{1}        = rho;
population.Kparams{1}    = nan;
population.tolSMC(1)     = max(rho);
SMCparams                = SMCparams(xorder,:);
if strcmp(kernel,'gamma_MGK')
    SMCparamsLog         = SMCparamsLog(xorder,:);
end

%% II: SMC step
% Main loop
tstartSMC = tic;
while (currTol~=tolFinal)&&(population.acceptRate(end)>pmin_acc)
    
    % Advance one iteration
    popIdx = popIdx + 1;
    
    % Define tolerance
    currTol = max([prctile(population.rho{popIdx-1},alpha),tolFinal])                              % This step ensures that the tolerance is decreasing and not larger than the target tolerance
        
    % Drop particles below the tolerance
    population.tolSMC(popIdx)     = currTol;
    dropedParticles               = find(rho>currTol);
    rho(dropedParticles)          = [];
    weights(dropedParticles)      = [];
    SMCparams(dropedParticles,:)  = [];
	HMM_params(dropedParticles)   = [];
    weights                       = weights/sum(weights);
    if strcmp(kernel,'gamma_MGK')
       SMCparamsLog(dropedParticles,:) = [];
    end
    
    % Start of next sampling step
    Nreplenish       = N-numel(weights);                                                            % Particles to replenish
    Nalive           = N-Nreplenish;
    weightsReplenish = zeros(Nreplenish,1);
    SMCreplenish     = zeros(Nreplenish,indexes(end));
    SMCreplenishLog  = zeros(Nreplenish,indexes(end));
	HMMreplenish{Nreplenish,1} = [];
    rhoReplenish     = zeros(Nreplenish,1);
    
    % Compute kernel parameters for the joint
    if strcmp(kernel,'mc_DK')                                                                       % normalize full vector of parameters for multi-component perturbation
        SMCparams_original = SMCparams;
        SMCparams = SMCparams./repmat(sum(SMCparams,2),1,size(SMCparams,2));
        
    elseif strcmp(kernel,'gamma_MGK')                                                               % work on the log-transformed parameters
        SMCparams_original = SMCparams;
        SMCparams = SMCparamsLog;
        
    elseif strcmp(kernel,'MGK')                                                                     % work only with free parameters
        SMCparams_original = SMCparams;
        SMCparams = mappingFxn(SMCparams);
        
    elseif strcmp(kernel,'ilrt_MGK')                                                                % work only with free parameters on the transformed space
        SMCparams_original = SMCparams;
        SMCparams = log(SMCparams)*fullMapping;
        
    elseif strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')
        SMCparams_original = SMCparams;
        SMCtemp = [];
        for ix = 1:numDirParams
            SMCtemp = [SMCtemp,mappingFxn_beta(SMCparams(:,indexes(ix,1):indexes(ix,2)))];
        end
        SMCparams = SMCtemp;
        if strcmp(kernel,'beta_logitMGK')
            SMCparams = logitMapping(SMCparams);
        end
        clearvars('SMCtemp');
    end
    
    % Compute kernel parameters
    Kparam = getGeneralKernelParam(SMCparams,kernel,method,weights,alpha,indexes,numDirParams);
    population.Kparams{popIdx} = Kparam;
    
    % Sampling loop
    counter    = 1;
    totalCount = 0;
    
    while (counter<=Nreplenish)
        
        % Update total counter
        totalCount = totalCount+1;
        
        % Resample from previous population
        prevIdx      = sum(rand(1)>=cumsum([0,weights']));
        prevParticle = SMCparams(prevIdx,:);
        
        % Perform perturbation according to the kernel
        if strcmp(kernel,'cw_DK')                                                           % perform DK perturbation on each Dir param
            newProbTheta{numDirParams,1} = [];
            thetaProposal = [];
            for ix = 1:numDirParams                
                probProposal     = randg(Kparam(ix)*prevParticle(indexes(ix,1):indexes(ix,2))+1);
                probProposal     = probProposal/sum(probProposal);
                thetaProposal    = [thetaProposal,probProposal];
                newProbTheta{ix} = probProposal;
            end            

        elseif strcmp(kernel,'mc_DK')                                                                       % perform DK perturbation on the Dir composition                       
            newParticle = randg(Kparam*prevParticle+1);
            newParticle = newParticle/sum(newParticle);
            newProbTheta{numDirParams,1} = [];
            thetaProposal = [];
            for ix = 1:numDirParams 
                newProbTheta{ix} = newParticle(indexes(ix,1):indexes(ix,2));
                newProbTheta{ix} = newProbTheta{ix}/sum(newProbTheta{ix});
                thetaProposal    = [thetaProposal,newProbTheta{ix}];
            end
            
        elseif strcmp(kernel,'MGK')                                                                         % perform MGK on the Dir composition
            newParticle = mvnrnd(prevParticle,Kparam);
            newParticle = invMapping(newParticle);
            newProbTheta{numDirParams,1} = [];
            thetaProposal = [];
            for ix = 1:numDirParams
                newProbTheta{ix} = newParticle(indexes(ix,1):indexes(ix,2));
                thetaProposal    = [thetaProposal,newProbTheta{ix}];
            end

        elseif strcmp(kernel,'ilrt_MGK')||strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')         % perform MGK on the transformed space
            newParticle = mvnrnd(prevParticle,Kparam);
            newProbTheta{numDirParams,1} = [];
            thetaProposal = [];
            for ix = 1:numDirParams
                if ~strcmp(kernel,'beta_logitMGK')
                    newProbTheta{ix} = invMapping{ix}(newParticle(transIdxs(ix,1):transIdxs(ix,2)));
                else
                    ybeta = invlogitMapping(newParticle(transIdxs(ix,1):transIdxs(ix,2)));
                    newProbTheta{ix} = invMapping{ix}(ybeta);
                end
                thetaProposal = [thetaProposal,newProbTheta{ix}];
            end
            
        elseif strcmp(kernel,'gamma_MGK')
            newParticle      = mvnrnd(prevParticle,Kparam);
            thetaProposalLog = newParticle;
            newProbTheta{numDirParams,1} = [];
            thetaProposal = [];            
            for ix = 1:numDirParams                
                newProbTheta{ix} = invMapping(newParticle(indexes(ix,1):indexes(ix,2)));                                
                thetaProposal    = [thetaProposal,newProbTheta{ix}];
            end
        end
        
        % Checkpoint for feasibility of the proposal
        if any(thetaProposal<0); continue; end; 
        simNewProbTheta = reshapeProbs(newProbTheta,priors,numDirParams,'simulation');
        
        % Accept/Reject particle
        simulated = getSummaryStatistics(simulateBayesNet(simNewProbTheta,model.n));
        rhoNew = sqrt(sum((simulated(:)-model.Sdata(:)).^2)/sum(model.Sdata(:)~=0));                                                % Calculate discrepancy                                                                                       
        
        % Check accuracy of the new particle
        if (rhoNew<=currTol)
            SMCreplenish(counter,:) = thetaProposal;            
            rhoReplenish(counter)   = rhoNew;
			HMMreplenish{counter,1} = simNewProbTheta;
            if strcmp(kernel,'gamma_MGK')
                SMCreplenishLog(counter,:) = thetaProposalLog;
            end
            
            % Compute (log)weight for the accepted particle
            if strcmp(kernel,'cw_DK')||strcmp(kernel,'mc_DK')||strcmp(kernel,'MGK')
                probPrior = 0;
                for ix = 1:numDirParams
                    probPrior = probPrior + probDirichlet(samplingPriors{ix},newProbTheta{ix},1);                                                  % These are log-probabilities
                end
                
                if strcmp(kernel,'cw_DK')
                    particleWeight = log(weights);
                    for ix = 1:numDirParams
                        particleWeight = particleWeight + probDirichlet(Kparam(ix)*SMCparams(:,indexes(ix,1):indexes(ix,2))+1,newProbTheta{ix},1);
                    end                    
                    particleWeight = logsumexp(particleWeight);
                    
                elseif strcmp(kernel,'mc_DK')
                    particleWeight = logsumexp(log(weights) + probDirichlet(Kparam*SMCparams+1,newParticle,1));
                    
                elseif strcmp(kernel,'MGK')
                    particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,mappingFxn(newParticle),Kparam)));
                end
                
            elseif strcmp(kernel,'ilrt_MGK')
                probPrior = 0;
                for ix = 1:numDirParams
                    probPrior  = probPrior + log_pz{ix}(mappingFxn{ix}(newProbTheta{ix}),samplingPriors{ix});                                             % These are log-probabilities
                end 
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,newParticle,Kparam)));
                
            elseif strcmp(kernel,'beta_MGK')
                probPrior = 0;
                for ix = 1:numDirParams                    
                    newProposal = mappingFxn_beta(newProbTheta{ix});                    
                    for jx = 1:numel(newProposal)
                        probPrior = probPrior + probDirichlet(newPriors{ix}(jx,:),[newProposal(jx),1-newProposal(jx)],1);   % Note that the Beta distribution is the simplest Dirichlet distribution
                    end
                end
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,newParticle,Kparam)));
                
            elseif strcmp(kernel,'beta_logitMGK')
                probPrior   = 0;
                newProposal = [];
                for ix = 1:numDirParams
                    temp        = logitMapping(mappingFxn_beta(newProbTheta{ix}));
                    probPrior   = probPrior + sum(log_pz(temp,newPriors{ix}));
                    newProposal = [newProposal,temp];
                end
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,newProposal,Kparam)));
                
            elseif strcmp(kernel,'gamma_MGK')                
                probPrior   = 0;
                newProposal = [];
                for ix = 1:numDirParams
                    probPrior   = probPrior + sum(log_pz(gammaMapping(newProbTheta{ix}),newPriors{ix}));
                end                                
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,thetaProposalLog,Kparam)));
            end
            
            % Compute un-normalized log-weight of replenished particle
            weightsReplenish(counter) = probPrior - particleWeight;
            
            % Update (accepted) particle counter
            counter = counter+1;
        end
    end
    
    % Re-normalize importance weights for robust calculations
    weightsReplenish = weightsReplenish - max(weightsReplenish);
    weightsReplenish = exp(weightsReplenish);
    
    % Normalize weights for the all the particles
    weightsReplenish = weightsReplenish/sum(weightsReplenish);
    weights          = [Nalive*weights;Nreplenish*weightsReplenish]/N;
    rho              = [rho;rhoReplenish];
	HMM_params       = [HMM_params;HMMreplenish];
    
    % Append new parameters according to their size
    if strcmp(kernel,'mc_DK')||strcmp(kernel,'MGK')||strcmp(kernel,'ilrt_MGK')||strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')
        SMCparams = [SMCparams_original;SMCreplenish];
    elseif strcmp(kernel,'gamma_MGK')        
        SMCparamsLog = [SMCparamsLog;SMCreplenishLog];
        SMCparams    = [SMCparams_original;SMCreplenish];
    else
        SMCparams = [SMCparams;SMCreplenish];
    end
    
    % Sort particles by tolerance score
    [rho,xorder] = sort(rho,'ascend');    
    SMCparams    = SMCparams(xorder,:);
	HMM_params   = HMM_params(xorder,:);
    weights      = weights(xorder);
    if strcmp(kernel,'gamma_MGK')
        SMCparamsLog = SMCparamsLog(xorder,:);
    end
    
    % Make assignation
    population.weights{popIdx}    = weights;
    population.acceptRate(popIdx) = (counter-1)/totalCount;
    population.ESS(popIdx)        = 1/sum(weights.^2);
    population.rho{popIdx}        = rho;
end

% Collect final samples
population.time   = toc(tstartSMC);
population.method = method;
population.kernel = kernel;

if (currTol==tolFinal)
    population.SMCparams = SMCparams;
	population.HMMparams = HMM_params;
else
    disp(['Acceptance rate below minimum, kernel: ',kernel,'.']);
end
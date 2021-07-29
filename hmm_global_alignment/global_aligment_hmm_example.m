function population = global_aligment_hmm_example(model,priors,N,alpha,tolStart,tolFinal,kernel,method)
% This scripts runs the single-reaction model
% Inputs:
%
% Outputs: population (structure array)
%--------------------- Pedro Saa UC 2021 ----------------------------------
% Initializing sampling parameters
pmin_acc   = 1e-5;                                                               % min. acceptance probability
currTol    = tolStart;                                                           % initial tolerance
SMCparams  = zeros(N,numel(priors{1})+numel(priors{2}));
rho        = zeros(N,1);
ESS        = N;
weights    = ones(N,1)/N;
numPriors  = numel(priors);
HMM_params = zeros(N,numel(priors{1})+numel(priors{2}));

% Define auxiliary functions
% x(1) = gamma , x(2) = delta , x(3) = nu , x(4) = epsil , x(5) = tau
TRANS    = @(x) [x(4),x(3),0,x(5);x(2),x(1),x(2),x(5);0,x(3),x(4),x(5);zeros(1,4)];
hmm_pars = @(theta1,theta2) [theta1(:,1).*(1-theta2(:,3)),theta1(:,2).*(1-theta2(:,3))/2,theta2];

% Initalize special transformation features for particular kernels
if strcmp(kernel,'gamma_MGK')
    
    % Initialize array with log( SMCparams )
    SMCparamsLog = zeros(N,numel(priors{1})+numel(priors{2}));
    
    % Construct prior distributions    
    for ix = 1:numPriors
        newPriors{ix} = priors{ix};        
        
        % Write mapping and inverse maping functions
        if (ix==1)
            invMapping   = @(y) exp(y)/sum(exp(y));
            gammaMapping = @(x) log(x);
            
            % Formulate transformed log-prior density
            %               k*z  - exp(z) - log( gamma( k ) )
            log_pz = @(z,k) k(:).*z(:) - exp(z(:)) - gammaln(k(:));
            
            % Initialize parameters for the Markov Chain (sampling directly from the prior is not possible in this case)
            burnin   = 2e3;
            thinning = 5e1;
        end
        
        % Generate a suitable covariance matrix by sampling from the
        % original prior, mapping the variables onto the transformed space,
        % and then computing the sample covariance matrix
        x          = randg(newPriors{ix}(ones(1e4,1),:));
        y          = gammaMapping(x);
        Q{ix}      = diag(cov(y));                                         % transtion kernel N(0,sigma) (symmetric)
        z_curr{ix} = mean(y);
    end
    
elseif strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')
    
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
    fullMapping           = [];
    U_full{numPriors}     = [];
    mappingFxn{numPriors} = [];
    invMapping{numPriors} = [];
    Q{numPriors}      = [];
    z_curr{numPriors} = [];
    detJ{numPriors}   = [];
    log_pz{numPriors} = [];
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
        
        % Formulate transformed log-prior density (we map z onto the full space)
        %               log( exp(             sum( )       )         )                      +  log( det_term ) -     log ( betaFxn )
        log_pz{ix} = @(z,alpha) (alpha-1)*(U_full{ix}*(log(invMapping{ix}(z))*U_full{ix})') + log(detJ{ix}(z)) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));
        
        % Initialize state of the Markov Chain
        z_curr{ix} = mappingFxn{ix}(priors{ix}/sum(priors{ix}));
        
        % Generate a suitable covariance matrix by sampling from the
        % original prior, mapping the variables onto the transformed space,
        % and then computing the sample covariance matrix
        x = randg(priors{ix}(ones(1e4,1),:));
        x = x./repmat(sum(x,2),1,numel(priors{ix}));
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
            z1_cand = mvnrnd(z_curr{1},Q{1});
            z2_cand = mvnrnd(z_curr{2},Q{2});
            
            % Compute acceptance prob
            accProb = exp( log_pz{1}(z1_cand,priors{1}) + log_pz{2}(z2_cand,priors{2}) - ...
                (log_pz{1}(z_curr{1},priors{1}) + log_pz{2}(z_curr{2},priors{2})) );
            
            % Accept proposal
            if (rand(1) <= accProb)
                z_curr{1} = z1_cand;
                z_curr{2} = z2_cand;
            end
        end
        
        % Map proposals to the simplex
        theta1 = invMapping{1}(z_curr{1});
        theta2 = invMapping{2}(z_curr{2});
        
        % Sample from beta/logit or log-gamma prior if required    
    elseif strcmp(kernel,'beta_logitMGK')||strcmp(kernel,'gamma_MGK')
        
        % Start MCMC from the previous point and check if burn-in is
        % required
        if (burnin>0)
            numSteps = burnin;
            burnin   = 0;                              % No more burn-in
            if strcmp(kernel,'beta_logitMGK')
                newPriors_pz = cell2mat(newPriors(:));
            else
                newPriors_pz = cell2mat(newPriors);
            end
        else
            numSteps = thinning;
        end
        
        % Run Metropolis-Hastings
        for ix = 1:numSteps
            z1_cand = mvnrnd(z_curr{1},diag(Q{1}));
            z2_cand = mvnrnd(z_curr{2},diag(Q{2}));
            
            % Compute acceptance prob
            accProb = exp( sum(log_pz([z1_cand,z2_cand],newPriors_pz)) - sum(log_pz([z_curr{1},z_curr{2}],newPriors_pz)) );
            
            % Accept proposal
            if (rand(1) <= accProb)
                z_curr{1} = z1_cand;
                z_curr{2} = z2_cand;
            end
        end
        
        % Map proposals to the simplex
        if strcmp(kernel,'beta_logitMGK')
            ybeta{1} = invlogitMapping(z_curr{1});
            ybeta{2} = invlogitMapping(z_curr{2});
            theta1   = invMapping{1}(ybeta{1});
            theta2   = invMapping{2}(ybeta{2});
        elseif strcmp(kernel,'gamma_MGK')
            thetaLog1 = z_curr{1};
            thetaLog2 = z_curr{2};
            theta1    = invMapping(thetaLog1);
            theta2    = invMapping(thetaLog2);
        end
        
        % Sample from beta prior if required
    elseif strcmp(kernel,'beta_MGK')
        
        % Draw from (D-1) independent beta distributions for each prior
        % Sample reversibilities
        theta1 = randg(newPriors{1});
        theta1 = (theta1(:,1)./sum(theta1,2))';
        
        % Sample enzyme abundances
        theta2 = randg(newPriors{2});
        theta2 = (theta2(:,1)./sum(theta2,2))';
        
        % Map proposals back to the simplex
        theta1 = invMapping{1}(theta1);
        theta2 = invMapping{2}(theta2);

        % Sample from the original prior for remaining cases
    else
        
        % Sample reversibilities
        theta1 = randg(priors{1});
        theta1 = theta1/sum(theta1);
        
        % Sample enzyme abundances
        theta2 = randg(priors{2});
        theta2 = theta2/sum(theta2);
    end
    
    % Calculate HMM parameters
    hmm_theta = hmm_pars(theta1,theta2);
    trans     = TRANS(hmm_theta);
    simulated = hist(hmmgenerate(model.n,model.TRANS(trans),model.EMIS));

    % Accept/Reject particle
    rhoNew = sum(abs(simulated-model.xdata));
    if (rhoNew<=currTol)
        rho(counter)          = rhoNew;
        SMCparams(counter,:)  = [theta1,theta2];
		HMM_params(counter,:) = hmm_theta;
        if strcmp(kernel,'gamma_MGK')
            SMCparamsLog(counter,:) = [thetaLog1,thetaLog2];
        end        
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
SMCparamsLog             = SMCparamsLog(xorder,:);
end

%% II: SMC step
% Main loop
tstartSMC = tic;
while (currTol~=tolFinal)&&(population.acceptRate(end)>pmin_acc)
    
    % Advance one iteration
    popIdx = popIdx + 1;
    
    % Define tolerance
    currTol = max([min([currTol-1,prctile(population.rho{popIdx-1},alpha)]),tolFinal]);             % This step ensures that the tolerance is decreasing and not larger than the target tolerance
    
    % Drop particles below the tolerance
    population.tolSMC(popIdx)     = currTol;
    dropedParticles               = find(rho>currTol);
    rho(dropedParticles)          = [];
    weights(dropedParticles)      = [];
    SMCparams(dropedParticles,:)  = [];
	HMM_params(dropedParticles,:) = [];
    weights                       = weights/sum(weights);
    if strcmp(kernel,'gamma_MGK')
       SMCparamsLog(dropedParticles,:) = [];
    end
    
    % Start of next sampling step
    Nreplenish       = N-numel(weights);                                                            % Particles to replenish
    Nalive           = N-Nreplenish;
    weightsReplenish = zeros(Nreplenish,1);
    SMCreplenish     = zeros(Nreplenish,numel(priors{1})+numel(priors{2}));
	HMMreplenish     = zeros(Nreplenish,numel(priors{1})+numel(priors{2}));
    rhoReplenish     = zeros(Nreplenish,1);
    
    % Compute kernel parameters for the joint
    if strcmp(kernel,'mc_DK')                                                                       % normalize full vector of parameters for multi-component perturbation
        SMCparams = SMCparams./repmat(sum(SMCparams,2),1,size(SMCparams,2));
        
    elseif strcmp(kernel,'gamma_MGK')                                                               % work on the log-parameters
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
        prevIdx = 1;
        for ix = 1:numPriors
            currIdx = prevIdx + numel(priors{ix}) - 1;
            if (ix==1)
                SMCtemp = mappingFxn_beta(SMCparams(:,prevIdx:currIdx));
            else
                SMCtemp = [SMCtemp,mappingFxn_beta(SMCparams(:,prevIdx:currIdx))];
            end
            prevIdx = currIdx + 1;
        end
        SMCparams = SMCtemp;
        if strcmp(kernel,'beta_logitMGK')
            SMCparams = logitMapping(SMCparams);            
        end
        clearvars('SMCtemp');
    end
    
    % Compute kernel parameters
    Kparam = getGeneralKernelParam(SMCparams,priors,kernel,method,weights,alpha);
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
            prevIdx = 1;
            for ix = 1:numPriors
                currIdx     = prevIdx + numel(priors{ix}) - 1;
                newParticle = randg(Kparam(ix)*prevParticle(prevIdx:currIdx)+1);
                if (ix==1)
                    theta1_new = newParticle/sum(newParticle);
                else
                    theta2_new = newParticle/sum(newParticle);
                end
                prevIdx = currIdx + 1;
            end
            newParticle = [theta1_new,theta2_new];
            
        elseif strcmp(kernel,'mc_DK')                                                       % perform DK perturbation on the Dir composition
            newParticle = randg(Kparam*prevParticle+1);
            newParticle = newParticle/sum(newParticle);
            prevIdx     = 1;
            for ix = 1:numPriors
                currIdx = prevIdx + numel(priors{ix}) - 1;
                if (ix==1)
                    theta1_new = newParticle(prevIdx:currIdx)/sum(newParticle(prevIdx:currIdx));
                else
                    theta2_new = newParticle(prevIdx:currIdx)/sum(newParticle(prevIdx:currIdx));
                end
                prevIdx = currIdx + 1;
            end
            
        elseif strcmp(kernel,'MGK')                                                         % perform MGK on the Dir composition
            newParticle = mvnrnd(prevParticle,Kparam);
            newParticle = invMapping(newParticle);
            prevIdx     = 1;
            for ix = 1:numPriors
                currIdx = prevIdx + numel(priors{ix}) - 1;
                if (ix==1)
                    theta1_new = newParticle(prevIdx:currIdx);
                else
                    theta2_new = newParticle(prevIdx:currIdx);
                end
                prevIdx = currIdx + 1;
            end
            
        elseif strcmp(kernel,'ilrt_MGK')||strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')   % perform MGK on the transformed space
            newParticle = mvnrnd(prevParticle,Kparam);
            prevIdx     = 1;
            for ix = 1:numPriors
                currIdx = prevIdx + numel(priors{ix}) - 2;
                if ~strcmp(kernel,'beta_logitMGK')
                    if (ix==1)
                        theta1_new = invMapping{ix}(newParticle(prevIdx:currIdx));
                    else
                        theta2_new = invMapping{ix}(newParticle(prevIdx:currIdx));
                    end
                else                    
                    ybeta = invlogitMapping(newParticle(prevIdx:currIdx));                                       
                    if (ix==1)
                        theta1_new = invMapping{ix}(ybeta);
                    else
                        theta2_new = invMapping{ix}(ybeta);
                    end                
                end
                prevIdx = currIdx + 1;
            end
            newParticle = [theta1_new,theta2_new];
            
        elseif strcmp(kernel,'gamma_MGK')
            newParticle = mvnrnd(prevParticle,Kparam);
            prevIdx     = 1;
            for ix = 1:numPriors
                currIdx = prevIdx + numel(priors{ix}) - 1;
                newProposal{ix} = newParticle(prevIdx:currIdx);
                if (ix==1)
                    theta1_new = invMapping(newProposal{ix});
                else
                    theta2_new = invMapping(newProposal{ix});
                end
                prevIdx = currIdx + 1;
            end
            newParticle = [theta1_new,theta2_new];
        end
        
        % Checkpoint for feasibility of the proposal
        if any(newParticle<0); continue; end;
        
        % Calculate HMM parameters
        hmm_theta = hmm_pars(theta1_new,theta2_new);
        trans     = TRANS(hmm_theta);
        simulated = hist(hmmgenerate(model.n,model.TRANS(trans),model.EMIS));
        
        % Accept/Reject particle
        rhoNew = sum(abs(simulated-model.xdata));                                                                                   % Calculate discrepancy                                                                                       
        
        % Check accuracy of the new particle
        if (rhoNew<=currTol)
            SMCreplenish(counter,:) = newParticle;
            rhoReplenish(counter)   = rhoNew;
			HMMreplenish(counter,:) = hmm_theta;
            if strcmp(kernel,'gamma_MGK')
                SMCreplenishLog(counter,:) = [newProposal{1},newProposal{2}];
            end
            
            % Compute (log)weight for the accepted particle
            if strcmp(kernel,'cw_DK')||strcmp(kernel,'mc_DK')||strcmp(kernel,'MGK')
                probPrior = probDirichlet(priors{1},theta1_new,1) + probDirichlet(priors{2},theta2_new,1);                          % These are log-probabilities
                
                if strcmp(kernel,'cw_DK')
                    particleWeight = log(weights);
                    prevIdx        = 1;
                    for ix = 1:numPriors
                        currIdx        = prevIdx + numel(priors{ix}) - 1;
                        particleWeight = particleWeight + probDirichlet(Kparam(ix)*SMCparams(:,prevIdx:currIdx)+1,newParticle(prevIdx:currIdx),1);
                        prevIdx        = currIdx + 1;
                    end
                    particleWeight = logsumexp(particleWeight);
                    
                elseif strcmp(kernel,'mc_DK')
                    particleWeight = logsumexp(log(weights) + probDirichlet(Kparam*SMCparams+1,newParticle,1));
                    
                elseif strcmp(kernel,'MGK')
                    particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,mappingFxn(newParticle),Kparam)));
                end
                
            elseif strcmp(kernel,'ilrt_MGK')
                probPrior      = log_pz{1}(mappingFxn{1}(theta1_new),priors{1}) + log_pz{2}(mappingFxn{2}(theta2_new),priors{2});   % These are log-probabilities
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,log(newParticle)*fullMapping,Kparam)));
                
            elseif strcmp(kernel,'beta_MGK')
                newProposal{numPriors} = [];
                probPrior = 0;
                for ix = 1:numPriors
                    if (ix==1)
                        newProposal{ix} = mappingFxn_beta(theta1_new);
                    else
                        newProposal{ix} = mappingFxn_beta(theta2_new);
                    end
                    for jx = 1:numel(newProposal{ix})
                        probPrior = probPrior + probDirichlet(newPriors{ix}(jx,:),[newProposal{ix}(jx),1-newProposal{ix}(jx)],1);   % Note that the Beta distribution is the simplest Dirichlet distribution
                    end
                end
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,[newProposal{1},newProposal{2}],Kparam)));
                
            elseif strcmp(kernel,'beta_logitMGK')
                for ix = 1:numPriors
                    if (ix==1)
                        newProposal{ix} = logitMapping(mappingFxn_beta(theta1_new));
                    else
                        newProposal{ix} = logitMapping(mappingFxn_beta(theta2_new));
                    end
                end
                probPrior = sum(log_pz([newProposal{1},newProposal{2}],newPriors_pz));
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,[newProposal{1},newProposal{2}],Kparam)));
                
            elseif strcmp(kernel,'gamma_MGK')          
                probPrior = sum(log_pz([newProposal{1},newProposal{2}],newPriors_pz));
                particleWeight = logsumexp(log(weights) + log(mvnpdf(SMCparams,[newProposal{1},newProposal{2}],Kparam)));
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
    if strcmp(kernel,'MGK')||strcmp(kernel,'ilrt_MGK')||strcmp(kernel,'beta_MGK')||strcmp(kernel,'beta_logitMGK')
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
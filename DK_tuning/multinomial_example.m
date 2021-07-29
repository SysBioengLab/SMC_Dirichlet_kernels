function population = multinomial_example(xdata,prior,N,alpha,tolStart,tolFinal,method)
% This scripts runs the multinomial example model
% Inputs:   % (1) MLE; (2) Total Var; (3) Gen Var
%
% Outputs: population (structure array) 
% --------------------- Pedro Saa UC 2021 ---------------------------------
%% Initialize parameters
clc
rng('default');
weights  = ones(N,1)/N;
thetaSMC = zeros(N,size(prior,2));
rho      = zeros(N,1);
currTol  = tolStart;
ESS      = N;

%% I: Prior sampling
tstartSMC  = tic;
counter    = 1;
totalCount = 1;
popIdx     = 1;
while (counter<=N)    
	
    % Sampling phase
    theta = randg(prior);
    theta = theta/sum(theta);
    
    % Simulate distribution
    xsim   = mnrnd(sum(xdata),theta);
    rhoNew = sum(abs(xsim-xdata));
    if (rhoNew<=currTol)
        thetaSMC(counter,:) = theta;
        rho(counter)        = rhoNew;
        counter             = counter+1;
    end
    
    % Update total counter
	totalCount = totalCount+1;
end

% Sort particles by tolerance score
[rho,xorder] = sort(rho,'ascend');

% Initialize population structure
thetaSMC                 = thetaSMC(xorder,:);
population.weights{1}    = weights(xorder);
population.acceptRate(1) = (counter-1)/(totalCount-1);
population.ESS(1)        = ESS;
population.rho{1}        = rho;
population.Kparams(1)    = nan;
population.tolSMC(1)     = max(rho);

%% II: SMC step
% Main loop
while (currTol~=tolFinal)
    
	% Advance one iteration
    popIdx = popIdx + 1;
    		
    % Define tolerance
    currTol = max([min([currTol-1,prctile(population.rho{popIdx-1},alpha)]),tolFinal]);           % This step ensures that the tolerance schedule is strictly decreasing
    
    % Drop particles below the tolerance
    population.tolSMC(popIdx)   = currTol;
    dropedParticles             = find(rho>currTol);
    rho(dropedParticles)        = [];
    weights(dropedParticles)    = [];
    thetaSMC(dropedParticles,:) = [];
    weights                     = weights/sum(weights);   
    
    % Start of next sampling step
    Nreplenish       = N-numel(weights);                                     % Particles to replenish
    Nalive           = N-Nreplenish;
    weightsReplenish = zeros(Nreplenish,1);
    thetaReplenish   = zeros(Nreplenish,size(prior,2));
    rhoReplenish     = zeros(Nreplenish,1);
    
    % Compute kernel parameters
    Kparam = getKernelParam(thetaSMC,method,weights,alpha);
    population.Kparams(popIdx) = Kparam;

    % Sequential loop
    counter    = 1;
    totalCount = 1;
    
    while (counter<=Nreplenish)
        
        % Resample from previous population
        prevIdx = sum(rand(1)>=cumsum([0,weights']));
        
        % Perform perturbation
        newParticle = randg(Kparam*thetaSMC(prevIdx,:)+1);
        newParticle = newParticle/sum(newParticle);
        
        % Simulate HMM
        xsim   = mnrnd(sum(xdata),newParticle);
        rhoNew = sum(abs(xsim-xdata));

        % Check accuracy of the new particle
        if (rhoNew<=currTol)         
            thetaReplenish(counter,:) = newParticle;
            rhoReplenish(counter)     = rhoNew;
            
            % Compute (log)weight for the accepted particle
            probPrior      = probDirichlet(prior,newParticle,1);
            particleWeight = logsumexp(log(weights) + probDirichlet(Kparam*thetaSMC+1,newParticle,1));
            weightsReplenish(counter) = probPrior - particleWeight;
                        
            % Update particle counter
            counter = counter+1;
        end
        
        % Update total counter
		totalCount = totalCount + 1;
    end
    
    % Re-normalize importance weights for robust calculations
    weightsReplenish = weightsReplenish - max(weightsReplenish);
    weightsReplenish = exp(weightsReplenish);

    % Normalize weights for the all the particles
    weightsReplenish = weightsReplenish/sum(weightsReplenish);
    weights          = [Nalive*weights;Nreplenish*weightsReplenish]/N;
    rho              = [rho;rhoReplenish];
    thetaSMC         = [thetaSMC;thetaReplenish];
    
    % Sort particles by tolerance score
    [rho,xorder] = sort(rho,'ascend');    
    thetaSMC     = thetaSMC(xorder,:);
    weights      = weights(xorder);
    
    % Make assignation    
    population.weights{popIdx}    = weights;    
    population.acceptRate(popIdx) = (counter-1)/(totalCount-1);
    population.ESS(popIdx)        = 1/sum(weights.^2);
    population.rho{popIdx}        = rho;
    if isnan(population.acceptRate(popIdx)); population.acceptRate(popIdx) = 1; end;            % This is just in case there is only one particle to replenish
end
population.time     = toc(tstartSMC);
population.method   = method;
population.thetaSMC = thetaSMC;
population.acceptRate(population.acceptRate==1) = [];
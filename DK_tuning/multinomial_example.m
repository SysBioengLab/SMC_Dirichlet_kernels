function population = multinomial_example(xdata,prior,N,alpha,tolStart,tolFinal,method)
% SMC-ABC with the Dirichlet kernel (DK) for a multinomial model, used by the
% DK-tuning study (Fig3.png, Supplementary Tables S1-S2).
%
% Inputs:   xdata     1 x D observed category counts
%           prior     1 x D Dirichlet prior hyper-parameters
%           N         number of particles
%           alpha     percentage of "alive" particles retained per iteration
%           tolStart  initial tolerance (use Inf)
%           tolFinal  target tolerance on the l1 count discrepancy
%           method    DK concentration tuning method, see getKernelParam:
%                     1 MLE concentration, 2 minimum KL divergence, 3 total variance
% Output:   population  scalar struct with fields weights, acceptRate, ESS, rho,
%                       Kparams, tolSMC (one entry per iteration), time,
%                       method and thetaSMC (the final N x D particles)
%
% The DK proposes inside the simplex by construction, so no proposal is ever
% rejected for falling outside it; rho is the l1 norm between simulated and
% observed counts, always an integer.
% --------------------- Pedro Saa UC 2026 ---------------------------------
%% Initialize parameters
rng('default');                                                            % every run starts from the same stream: the published results depend on it
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

% Sort particles by discrepancy
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
    population.tolSMC(popIdx)    = currTol;
    droppedParticles             = find(rho>currTol);
    rho(droppedParticles)        = [];
    weights(droppedParticles)    = [];
    thetaSMC(droppedParticles,:) = [];
    weights                      = weights/sum(weights);

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

        % Simulate counts from the multinomial
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

    % Normalize the weights over all the particles
    weightsReplenish = weightsReplenish/sum(weightsReplenish);
    weights          = [Nalive*weights;Nreplenish*weightsReplenish]/N;
    rho              = [rho;rhoReplenish];
    thetaSMC         = [thetaSMC;thetaReplenish];

    % Sort particles by discrepancy
    [rho,xorder] = sort(rho,'ascend');
    thetaSMC     = thetaSMC(xorder,:);
    weights      = weights(xorder);

    % Store this iteration
    population.weights{popIdx}    = weights;
    population.acceptRate(popIdx) = (counter-1)/(totalCount-1);
    population.ESS(popIdx)        = 1/sum(weights.^2);
    population.rho{popIdx}        = rho;
    if isnan(population.acceptRate(popIdx)); population.acceptRate(popIdx) = 1; end          % 0/0: no particle needed replenishing (Nreplenish = 0)
end
population.time     = toc(tstartSMC);
population.method   = method;
population.thetaSMC = thetaSMC;
% Drop the entries equal to 1: the prior stage, where every draw is accepted when
% tolStart = Inf, and any iteration with nothing to replenish.
population.acceptRate(population.acceptRate==1) = [];
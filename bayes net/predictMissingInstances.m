function predictMissingInstances(population)
%------------------------- Pedro Saa UC 2021 ------------------------------
clc

% Initialize model parameters
rng('default');                                                              % for reproducibility
load('probs.mat');

% Generate synthetic data
model.n        = 1e3;                                                        % dataset size
model.probCond = probs;                                                      % conditional probabilities
model.xdata    = simulateBayesNet(model.probCond,model.n);
indexes        = getParamIndex(probs);
numDirParams   = size(indexes,1);

% Corrupt data with blank data
noiseLevel = .3;
ix = sum(repmat(rand(numel(model.xdata),1),1,3) >= cumsum([zeros(numel(model.xdata),1),noiseLevel*ones(numel(model.xdata),1),(1-noiseLevel)*ones(numel(model.xdata),1)],2),2);
model.xmissing = model.xdata;
model.xmissing(ix==1) = 0;

% Predict instances
missingInstances = find(any(model.xmissing'==0)');
numPredictions   = numel(missingInstances);
predInstance     = zeros(numel(missingInstances),6);
obsInstance      = zeros(numel(missingInstances),6);
classObs{6}      = [];   
classPred{6}     = [];
instancePred     = 0;
for ix = 1:numPredictions
    ix
    % Report progress
    obsInstance(ix,:) = model.xdata(missingInstances(ix),:);
    testInstance       = model.xmissing(missingInstances(ix),:);
    predInstanceModel = zeros(numel(population.HMMparams),6);
    
    % Loop through all particles
    for jx = 1:numel(population.HMMparams)
        
        % Simulate from the current particle if possible, otherwise change
        % the format of the particle
        try
            predInstanceModel(jx,:) = posteriorProbBayesNet(population.HMMparams{jx}(:),testInstance);
        catch
            predInstanceModel(jx,:) = posteriorProbBayesNet(reshapeProbs(population.HMMparams{jx}(:),probs,numDirParams,'simulation'),testInstance);
        end        
    end

    % Predict instance and clean predInstanceModel array
    predInstance(ix,:) = population.weights{end}'*predInstanceModel;

    % Determine predicted class and compare with the observations
    predInstance(ix,:) = getDiscreteInstance(predInstance(ix,:));

    if (all(predInstance(ix,:)==obsInstance(ix,:)))
        instancePred = instancePred +1;
    end
    
    % Find the missing clases and build array
    idxClassMissing = find(model.xmissing(missingInstances(ix),:)==0);
    for kx = idxClassMissing
        classObs{kx}  = [classObs{kx},model.xdata(missingInstances(ix),kx)];
        classPred{kx} = [classPred{kx},predInstance(ix,kx)];
    end
end

% Save final results
instancePred = instancePred/numPredictions;
save('missingPrediction.mat','classObs','classPred','instancePred');

function reshapedProbs = reshapeProbs(condProbs,priors,numDirParams,mode)
% Convert conditional probabilities between the two layouts used by the Bayes net.
%
% Inputs:   condProbs     probabilities in the source layout
%           priors        node priors, which fix each node's table size
%           numDirParams  number of Dirichlet blocks (used by 'sampling' only)
%           mode          'simulation': one cell per BLOCK -> one table per NODE
%                         'sampling'  : one table per NODE -> one cell per BLOCK
% Output:   reshapedProbs the same probabilities in the target layout
% --------------------- Pedro Saa UC 2026 ---------------------------------
% Reshape cond probabilities for simulation (from sampling form)
if strcmp(mode,'simulation')
    reshapedProbs{numel(priors),1} = [];
    counter = 1;
    for ix  = 1:numel(priors)
        reshapedProbs{ix} = zeros(size(priors{ix}));
        for jx = 1:size(priors{ix},1)
            reshapedProbs{ix}(jx,:) = condProbs{counter};
            counter = counter + 1;
        end
    end

% Reshape cond probabilities for sampling (from simulation form)
elseif strcmp(mode,'sampling')
    reshapedProbs{numDirParams,1} = [];
    counter = 1;
    for ix  = 1:numel(condProbs)
        for jx = 1:size(condProbs{ix},1)
            reshapedProbs{counter} = condProbs{ix}(jx,:);
            counter = counter + 1;
        end
    end
end
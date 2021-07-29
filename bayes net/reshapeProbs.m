function reshapedProbs = reshapeProbs(condProbs,priors,numDirParams,mode)
% Reshape conditional probabilities into the proper form
% --------------------- Pedro Saa UC 2021 ---------------------------------
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
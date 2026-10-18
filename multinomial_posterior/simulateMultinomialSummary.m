function S = simulateMultinomialSummary(theta,model)
% Draw model.n multinomial trials with probabilities theta and return the
% category counts.
%
% Inputs:   theta  1 x D category probabilities (a point on the simplex)
%           model  struct with field n (number of trials)
% Output:   S      1 x D counts, summing to model.n
%
% The counts are sufficient, so at a tolerance of zero the ABC posterior is the
% exact posterior Dir(prior + xdata) and any discrepancy from it is a property of
% the sampler alone. Written out rather than through mnrnd (an M-file) because it
% is called once per proposed particle on the sampler's hot path.
% --------------------- Pedro Saa UC 2026 ---------------------------------
D = numel(theta);
c = cumsum(theta(:)');
idx = 1 + sum(rand(model.n,1) >= c,2);
idx = min(idx,D);                                                          % guard against round-off in the cumulative sum
S   = accumarray(idx,1,[D 1])';

function S = simulateCovidSummary(theta,off,n)
% Simulate a cohort from the COVID-19 screening Bayes net and summarise it in one
% pass, without materialising the cohort.
%
% Inputs:   theta  1 x 72 flat parameter vector (all Dirichlet blocks)
%           off    node offsets from getCovidOffsets
%           n      number of subjects to simulate
% Output:   S      1 x 31 summary statistics (layout in getCovidSummary)
%
% Fuses three steps that were ~82% of the sampler's runtime: building the CPD cell
% array, filling an n x 14 cohort, and rescanning it. theta is read through the
% offsets, so no CPD tables are built.
% --------------------- Pedro Saa UC 2026 ---------------------------------
st = off.start;
cl = off.cols;

% Root node: COVID-19 status. m1/m2 are the two strata used by every statistic.
covState = (rand(n,1) >= theta(st(1)));                                          % false -> state 1, true -> state 2
m1  = ~covState;                 n1 = max(sum(m1),1);                            % COVID-negative; guard against an empty stratum
m2  =  covState;                 n2 = max(sum(m2),1);                            % COVID-positive
covIdx = double(covState);                                                       % 0 for state 1, 1 for state 2

S    = zeros(1,31);
S(1) = sum(m2)/n;                                                           % prevalence, i.e. the positive class (see getCovidSummary)

% Self-reported symptoms, conditional on the COVID-19 status
for j = 1:7
    k = 1+j;
    p = reshape(theta(st(k) + covIdx*cl(k)),n,1);                           % P(symptom = state 1 | COV)
    s = (rand(n,1) < p);
    S(1+j) = sum(s & m1)/n1;
    S(8+j) = sum(s & m2)/n2;
end

% Hidden olfactory status (3 states) given COV. The index is n x 3 so the gather
% returns n x 3; an n x 1 index would return a ROW and broadcast to n x 3n.
base = st(9) + covIdx*cl(9);
P    = theta([base,base+1,base+2]);
kor  = 1 + sum(rand(n,1) >= cumsum(P,2),2);
kor  = min(max(kor,1),3);                                                   % guard against round-off in the cumulative sum

% Psychophysical odour outcomes, conditional on the hidden status
K = zeros(n,1);
for k = 1:6
    nd = 9+k;
    p  = reshape(theta(st(nd) + (kor-1)*cl(nd)),n,1);                       % P(odour = state 1 | sKor)
    o  = (rand(n,1) < p);
    S(15+k) = sum(o & m1)/n1;
    S(21+k) = sum(o & m2)/n2;
    K = K + o;
end

% Central moments of the odour count per stratum. sum(x)/numel(x) rather than
% mean(x): mean is an M-file and this runs once per proposed particle.
K1 = K(m1);   c1 = numel(K1);   d1 = K1 - sum(K1)/c1;
K2 = K(m2);   c2 = numel(K2);   d2 = K2 - sum(K2)/c2;
S(28) = sum(d1.^2)/c1;   S(29) = sum(d1.^3)/c1;
S(30) = sum(d2.^2)/c2;   S(31) = sum(d2.^3)/c2;

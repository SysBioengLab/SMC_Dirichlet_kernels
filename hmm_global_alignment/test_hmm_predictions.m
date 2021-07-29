function test_hmm_predictions
%--------------------- Pedro Saa UC 2021 ----------------------------------
clc,clear
addpath('results')

% Initialize model parameters
rng('default');                                     % for reproducibility

% load parameters
load('results_hmm_full3.mat');                      % pick solution from ilrt kernel
thetaSMC = population_MGK.HMMparams;
weights  = population_MGK.weights{end};

% HMM standar parameters
p0     = [.2,.5,.2,.1];                             % initial probability parameters
EMIS   = [.25,zeros(1,5);...
          0,.5,.3,.15,.05,0;...                                            % emission probability matrix
          zeros(1,5),.25;...
          zeros(1,6)];       

% x(1) = gamma , x(2) = delta , x(3) = nu , x(4) = epsil , x(5) = tau
transFxn = @(x) [x(4),x(3),0,x(5);x(2),x(1),x(2),x(5);0,x(3),x(4),x(5);zeros(1,4)];
     
% Build probability matrices and generate synthetic data
EMIS     = [zeros(1,size(EMIS,2));EMIS];

% Nomenclature: 1: x insertion, 2: exact match, 3: AT/CG mismatch, 4: GT/AC mismatch, 5: CT/AG mismatch, 6: y insertion                    
% Emission seqs:                x: TTACG  ,   y = TAG
seqs = [2,6,2,6,2;...           %  T-A-G      : possible y alignments
        6,2,2,6,2;...           %  -TA-G
        6,6,3,4,2;...           %  --TAG
        2,3,5,6,6;...           %  TAG--
        6,2,6,4,2;...           %  -T-AG
        2,6,6,4,2;...           %  T--AG
        6,2,2,3,6;...           %  -TAG-
        2,6,2,3,6;...           %  T-AG-
        2,3,6,3,6;...           %  TA-G-
        2,3,6,6,2];             %  TA--G

% Loop through all the particles and sequences and compute posterior
% probabilities
M = size(seqs,1);
N = size(thetaSMC,1);
pseqs = zeros(M,1);
pstates{M} = [];
for ix = 1:M
    pseq = 0;
    pstates{ix} = 0;
    for jx = 1:N
        trans = transFxn(thetaSMC(jx,:));
        TRANS = [0,p0;zeros(size(trans,1),1),trans];
        [pstate,logpseq] = hmmdecode(seqs(ix,:),TRANS,EMIS);
        pseq        = pseq + exp(logpseq)*weights(jx);
        pstates{ix} = pstates{ix} + pstate*weights(jx);      
    end    
    pseqs(ix) = pseq;
end

% Get posterior probability using true model
trans = transFxn([.5,.2,.8,.1,.1]);
TRANS = [0,p0;zeros(size(trans,1),1),trans];
logpseq = zeros(M,1);
for ix = 1:M
    [~,logpseq(ix)] = hmmdecode(seqs(ix,:),TRANS,EMIS);
end

% Plot results
barh([pseqs/sum(pseqs),exp(logpseq)/sum(exp(logpseq))])
axis([0,.3,.5,10.5])
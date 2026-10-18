function priors = initializePriors_final_sim(fullData)
% Dirichlet priors of the 15 nodes of the COVID-19 screening Bayes net.
%
% Input:    fullData  n x 14 cohort with states coded false = 1, true = 2:
%                     [COV, 7 self-reported symptoms, 6 odour outcomes]
% Output:   priors    1 x 15 cell; priors{k} is (parent states) x (own states),
%                     each row being one Dirichlet block, in node order
%
% The values are Dirichlet hyper-parameters, i.e. pseudo-counts, not normalised
% probabilities. Nodes whose parent is observed are data-informed: each entry is
% the number of subjects falling in that cell of the cohort. The hidden olfactory
% status sKor and the six odour outcomes it governs have flat priors.
%------------------------- Pedro Saa UC 2026 ------------------------------
covState = fullData(:,1);
priors   = cell(1,15);

% COV, the root: [negative, positive]
priors{1} = [sum(covState==1), sum(covState==2)];

% Self-reported symptoms, each with COV as parent. Rows are COV negative /
% positive and columns symptom absent / present. Cohort columns 2-8 hold, in
% order: cold, cough, fever, muscular pain, breathing difficulty, ageusia and
% anosmia.
for k = 2:8
    s = fullData(:,k);
    priors{k} = [sum(covState==1 & s==1), sum(covState==1 & s==2);
                 sum(covState==2 & s==1), sum(covState==2 & s==2)];
end

% sKor, the hidden objective olfactory status (anosmic, hyposmic, normosmic),
% with COV as parent: 2 x 3, uniform.
priors{9} = ones(2,3);

% The six odour outcomes oB, oC, oM, oO, oP and oV, each dichotomous with sKor as
% parent: 3 x 2, uniform.
for k = 10:15
    priors{k} = ones(3,2);
end
end

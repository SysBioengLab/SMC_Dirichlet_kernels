function [tau,Jbest] = youdenThreshold(p,y)
% Operating threshold maximising Youden's J = sensitivity + specificity - 1.
%
% Inputs:   p     n x 1 predicted probabilities of the positive class
%           y     n x 1 logical, true class
% Outputs:  tau   selected threshold
%           Jbest the attained J
%
% Candidates are the distinct predicted probabilities, so every attainable
% operating point is examined. Shared by classifier_performance (Table 1) and by
% make_figures_covid (Fig. 6F) so that both mark the same operating point.
% --------------------- Pedro Saa UC 2026 ---------------------------------
c     = unique(p(:));
c     = [c; max(c)+eps(max(c))];                                           % also allow "classify nobody positive"
nPos  = sum(y); nNeg = numel(y)-nPos;
Jbest = -Inf; tau = c(1);
for i = 1:numel(c)
    yh = (p >= c(i));
    J  = sum(yh & y)/nPos - sum(yh & ~y)/nNeg;                             % sens - (1-spec)
    if (J > Jbest); Jbest = J; tau = c(i); end
end
end

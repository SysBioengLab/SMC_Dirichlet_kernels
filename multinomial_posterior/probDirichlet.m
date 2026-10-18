function log_prob = probDirichlet(alpha,xdata,flag,weights)
% Dirichlet log-density at xdata (flag = 1) or weighted log-likelihood on the
% log-scale alpha' = log(alpha), for unconstrained optimization (flag = 2).
%--------------------- Pedro Saa UC 2026 ----------------------------------
if (flag==1)
    log_betaFxn = sum(gammaln(alpha),2) - gammaln(sum(alpha,2));
    log_prob    = sum((alpha-1).*log(xdata(ones(1,size(alpha,1)),:)),2) - log_betaFxn;

elseif (flag==2)
    alpha       = exp(alpha);
    log_betaFxn = sum(gammaln(alpha),2) - gammaln(sum(alpha,2));
    log_prob    = weights'*(sum((alpha(ones(size(xdata,1),1),:)-1).*log(xdata),2) - log_betaFxn);
end

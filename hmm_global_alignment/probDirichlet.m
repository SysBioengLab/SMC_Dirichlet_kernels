function log_prob = probDirichlet(alpha,xdata,flag,weights)
%--------------------------------------------------------------------------
% Compute log:probability (flag=1) or weighted log-likelihood (flag=2) of Dir distribution
%
% Inputs:       xdata (in rows), alpha , weights  , flag
%
% Outputs:      prob (f or log-likelihood f)
%--------------------- Pedro Saa UC 2021 ----------------------------------
% Return probability evaluated at xdata
if (flag==1)    
    log_betaFxn = sum(gammaln(alpha),2) - gammaln(sum(alpha,2));                    % Compute log:beta_fxn
    log_prob    = sum((alpha-1).*log(xdata(ones(1,size(alpha,1)),:)),2) - log_betaFxn;
    
% Return weighted log:likelihood for unconstrained optimization (alpha' ~ log(alpha)
elseif (flag==2)
    alpha       = exp(alpha);
    log_betaFxn = sum(gammaln(alpha),2) - gammaln(sum(alpha,2));                    % Compute log:beta_fxn               
    log_prob    = weights'*(sum((alpha(ones(size(xdata,1),1),:)-1).*log(xdata),2) - log_betaFxn);
end

function fobj = compute_kld(Kparam,thetaSMC_q,thetaSMC_t,weights_q,weights_t,Nalive_t,Nalive_q)
%--------------------------------------------------------------------------
% Formulate KL problem for Kparams estimation
%
% Inputs:
%
% Outputs:
%--------------------- Pedro Saa UC 2021 ----------------------------------
% Initialize parameters
Kparam      = exp(Kparam);
q           = zeros(Nalive_t,1);
alpha       = Kparam*thetaSMC_q + 1;
log_betaFxn = sum(gammaln(alpha),2) - gammaln(sum(alpha,2));
for ix = 1:Nalive_t
    theta    = thetaSMC_t(ix,:);
    log_prob = sum((alpha - 1).*log(theta(ones(1,Nalive_q),:)),2) - log_betaFxn;
    q(ix)    = logsumexp(log(weights_q) + log_prob);
end
fobj = weights_t'*q;
end
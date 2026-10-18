function K = getKernelParam(theta,method,weights,alpha)
% Concentration kappa_t of the Dirichlet kernel, from the weighted population.
%
% Inputs:   theta    N x D particles on the simplex, sorted by discrepancy
%           method   1 MLE concentration, 2 minimum KL divergence, 3 total variance
%           weights  N x 1 normalised importance weights
%           alpha    percentage of alive particles (used by method 2 only)
% Output:   K        concentration; the kernel is Dir(K*theta_i + 1)
%
% This is the version that produced the DK-tuning study. The other case-study
% folders use getDKparam, whose method 1 is solved by SQP with an analytic
% gradient and clamped to K >= 0; the two agree whenever the clamp is inactive.
% --------------------- Pedro Saa UC 2026 ---------------------------------
D = size(theta,2);
switch method
    case 1                      % Pseudo concentration parameter based on MLE estimate
        options  = optimset('TolX',1e-5,'TolFun',1e-6,'MaxIter',1e4,'MaxFunEvals',1e4,'Display','off');
        alpha0   = log(ones(1,D));
        alphaMLE = fminsearch(@(alpha) -probDirichlet(alpha,theta,2,weights),alpha0,options);
        K        = sum(exp(alphaMLE))-D;

    case 2                      % Parameter based on the KL-divergence between proposal and target
        options    = optimset('TolX',1e-5,'TolFun',1e-6,'MaxIter',1e4,'MaxFunEvals',1e4,'Display','off');
        Nalive_t   = ceil(numel(weights)*alpha/100);
        weights_t  = weights(1:Nalive_t);
        weights_t  = weights_t/sum(weights_t);
        weights_q  = weights(Nalive_t+1:end);
        weights_q  = weights_q/sum(weights_q);
        thetaSMC_t = theta(1:Nalive_t,:);
        thetaSMC_q = theta(Nalive_t+1:end,:);
        Nalive_q   = numel(weights)-Nalive_t;
        K          = fminsearch(@(alpha) -compute_kld(alpha,thetaSMC_q,thetaSMC_t,weights_q,weights_t,Nalive_t,Nalive_q),log(10),options);
        K          = exp(K)/D;

    case 3                      % Weighted total variance
        Var = zeros(D);
        for ix = 1:D
            for jx = 1:D
                if (ix~=jx)
                    Var(ix,jx) = var(log(theta(:,ix)./theta(:,jx)),weights);
                end
            end
        end
        K = 2*D/sum(Var(:));
end
function K = getDKparam(theta,method,weights,alpha)
% Concentration kappa_t of the Dirichlet kernel for one Dirichlet block.
%
% Inputs:   theta    N x D particles of the block, sorted by discrepancy
%           method   1 MLE concentration, 2 minimum KL divergence, 3 total variance
%           weights  N x 1 normalised importance weights
%           alpha    percentage of alive particles (used by method 2 only)
% Output:   K        concentration, >= 0; the kernel is Dir(K*theta_i + 1)
% --------------------- Pedro Saa UC 2026 ---------------------------------
% Options built once and cached: building an optimoptions object costs more than
% the solve it configures. Run "clear getDKparam" after editing them.
persistent fminconOpts
if isempty(fminconOpts)
    fminconOpts = optimoptions('fmincon','Algorithm','sqp','Display','off',...
                               'SpecifyObjectiveGradient',true,...
                               'OptimalityTolerance',1e-8,'StepTolerance',1e-10,...
                               'MaxIterations',200,'MaxFunctionEvaluations',1e3);
end

D = size(theta,2);
switch method
    case 1                      % Pseudo concentration parameter based on MLE estimate
        % SQP on the log-scale, alpha = exp(a), with the analytic gradient
        options  = fminconOpts;
        a0       = zeros(1,D);                                                          % a = log(alpha) = 0  <=>  alpha = 1
        lb       = log(1e-6)*ones(1,D);                                                 % generous box, keeps SQP well posed
        ub       = log(1e6)*ones(1,D);
        aMLE     = fmincon(@(a) negDirLogLik(a,theta,weights),a0,[],[],[],[],lb,ub,[],options);
        K        = sum(exp(aMLE))-D;

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

% The concentration must be non-negative: the DK hyper-parameter is
% alpha = K*theta + 1 and needs alpha > 0. Method 1 returns sum(alpha_MLE) - D,
% negative for a population more diffuse than uniform; K = 0 then gives the
% broadest admissible kernel, Dirichlet(1,...,1).
K = max(K,0);
end

%% ------------------------------------------------------------------------
function [f,g] = negDirLogLik(a,theta,weights)
% Negative weighted Dirichlet log-likelihood and gradient on the log-scale
% alpha = exp(a), i.e. -probDirichlet(a,theta,2,weights) plus its gradient.
%   L      = sum_j (alpha_j-1)*m_j - W*(sum_j gammaln(alpha_j) - gammaln(sum alpha))
%   dL/da_j= alpha_j*( m_j - W*(psi(alpha_j) - psi(sum_k alpha_k)) )
%   m_j = sum_i w_i*log(theta_ij),  W = sum_i w_i
alpha = exp(a);
W     = sum(weights);
m     = weights'*log(theta);                                               % 1 x D weighted mean log-composition
sA    = sum(alpha);

f = -( (alpha-1)*m' - W*(sum(gammaln(alpha)) - gammaln(sA)) );
if nargout > 1
    g = -( alpha.*(m - W*(psi(alpha) - psi(sA))) );
end
end
function K = getKernelParam(theta,method,weights,alpha)
%--------------------------------------------------------------------------
% Computes global parameter of the dirichlet kernel
%
% Inputs:       theta (Dirichlet variates)
%
% Outputs:      K (kernel param)
%--------------------- Pedro Saa UC 2021 ----------------------------------
D = size(theta,2);
switch method
    case 1						 % Pseudo concentration parameter based on MLE estimate
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
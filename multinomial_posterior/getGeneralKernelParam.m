function [Kparam,shrinkDelta] = getGeneralKernelParam(SMCparams,priors,kernel,method,weights,alpha)
% Kernel parameters: per-block Dirichlet parameters for the DK, or twice the
% weighted sample covariance (shrunk if near-singular) for the MGK family.
%
% Outputs:  Kparam      (double array)
%           shrinkDelta (0 if the raw covariance was used; otherwise the
%                        block-diagonal shrinkage intensity that was applied)
%--------------------- Pedro Saa UC 2026 ----------------------------------
shrinkDelta = 0;

% Check type of kernel
if strcmp(kernel,'cw_DK')
    numPriors = numel(priors);
    Kparam(numPriors) = 0;
    prevIdx = 1;
    for ix = 1:numPriors
        currIdx    = prevIdx + numel(priors{ix}) - 1;
        Kparam(ix) = getDKparam(SMCparams(:,prevIdx:currIdx),method,weights,alpha);     % compute DK parameters for each Dir param
        prevIdx    = currIdx + 1;
    end

else
    muX        = weights'*SMCparams;
    diffX      = SMCparams-muX(ones(size(SMCparams,1),1),:);
    diffW      = weights(:,ones(1,size(SMCparams,2))).*diffX;
    Kparam     = 2*diffW'*diffX/(1-sum(weights.^2));                                    % twice weighted sample covariance
    [Kparam,shrinkDelta] = shrinkToBlockDiag(Kparam,priors,size(SMCparams,2));
end
end

%% ------------------------------------------------------------------------
function [K,delta] = shrinkToBlockDiag(K,priors,dWork)
% Block-diagonal shrinkage K(delta) = (1-delta)*K + delta*blockdiag(K), applied
% only when the weighted sample covariance is near-singular: as the population
% concentrates, its distinct-ancestor count falls below the working dimension and
% K loses rank. The target is block-diagonal rather than spherical because the
% blocks are the Dirichlet parameters themselves and stay well estimated. delta
% interpolates between the full kernel (0) and a per-block one (1), damping only
% cross-block correlation. The returned K is used for BOTH proposal and weight.
condTrigger = 1e10;                                                                     % regard K as near-singular beyond this condition number
condTarget  = 1e8;                                                                      % shrink just enough to get below this

K     = (K+K')/2;                                                                       % remove round-off asymmetry
delta = 0;
if ~isNearSingular(K,condTrigger)
    return
end

% Block partition of the WORKING space. The transformed kernels (MGK, ILRT,
% beta-logit) drop one coordinate per block; gamma_MGK keeps them all.
dOrig = zeros(numel(priors),1);
for ix = 1:numel(priors)
    dOrig(ix) = numel(priors{ix});
end
if (dWork==sum(dOrig))
    blk = dOrig;
else
    blk = dOrig-1;
end
blk    = blk(blk>0);
stopI  = cumsum(blk);
startI = stopI-blk+1;

% Block-diagonal target: keep the within-block covariances, drop cross-block terms
B = zeros(size(K));
for b = 1:numel(blk)
    r = startI(b):stopI(b);
    B(r,r) = K(r,r);
end

% Smallest intensity on a fixed ladder that restores the target conditioning
for d = [0.01 0.02 0.05 0.10 0.20 0.35 0.50 0.75 1.00]
    Kd = (1-d)*K + d*B;
    Kd = (Kd+Kd')/2;
    if ~isNearSingular(Kd,condTarget)
        K     = Kd;
        delta = d;
        return
    end
end

% Even the block-diagonal target is degenerate (a block has collapsed to a
% point): fall back to it plus the minimal ridge that restores conditioning.
B     = (B+B')/2;
ev    = eig(B);
K     = B + max(max(ev)/condTarget,realmin)*eye(size(B));
delta = 1;
end

%% ------------------------------------------------------------------------
function tf = isNearSingular(K,kappa)
ev = eig((K+K')/2);
tf = (min(ev)<=0) || (min(ev) <= max(ev)/kappa);
end

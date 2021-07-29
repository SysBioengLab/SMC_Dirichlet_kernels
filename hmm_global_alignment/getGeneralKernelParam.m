function Kparam = getGeneralKernelParam(SMCparams,priors,kernel,method,weights,alpha)
% This scripts calculates the kernel parameters according to different
% methods
% Inputs:
%
% Outputs:  Kparam (cell array)
%--------------------- Pedro Saa UC 2021 ----------------------------------
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
    
elseif strcmp(kernel,'mc_DK')
    Kparam = getDKparam(SMCparams,method,weights,alpha);                                % compute DK parameters for the composition of Dir params

else
    muX        = weights'*SMCparams;
    diffX      = SMCparams-muX(ones(size(SMCparams,1),1),:);
    diffW      = weights(:,ones(1,size(SMCparams,2))).*diffX;
    Kparam     = 2*diffW'*diffX/(1-sum(weights.^2));                                    % twice weighted sample covariance
end
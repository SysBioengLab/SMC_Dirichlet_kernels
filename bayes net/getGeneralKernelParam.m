function Kparam = getGeneralKernelParam(SMCparams,kernel,method,weights,alpha,indexes,numDirParams)
% This scripts calculates the kernel parameters according to different
% methods
% Inputs:
%
% Outputs:  Kparam (double array)
% --------------------- Pedro Saa UC 2021 ---------------------------------
% Check type of kernel
if strcmp(kernel,'cw_DK')
    Kparam(numDirParams) = 0;
    for ix = 1:numDirParams        
        Kparam(ix) = getDKparam(SMCparams(:,indexes(ix,1):indexes(ix,2)),method,weights,alpha);     % compute DK parameters for each Dir param
    end
    
elseif strcmp(kernel,'mc_DK')
    Kparam = getDKparam(SMCparams,method,weights,alpha);                                % compute DK parameters for the composition of Dir params
    
else
    muX        = weights'*SMCparams;
    diffX      = SMCparams-muX(ones(size(SMCparams,1),1),:);
    diffW      = weights(:,ones(1,size(SMCparams,2))).*diffX;
    Kparam     = 2*diffW'*diffX/(1-sum(weights.^2));                                    % twice weighted sample covariance
end
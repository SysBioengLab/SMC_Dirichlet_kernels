function spec = getKernelSpec(kernelIdx)
% Map the kernel INDEX used throughout the paper onto everything a run needs.
%
% Input:    kernelIdx  1, 3, 4, 6 or 7 (see the table below); omit to get the
%                      full 1 x 5 array of specs
% Output:   spec       struct with fields idx, kernel, method, paperName,
%                      varName and resultFile
%
%   idx  code kernel     paper name   DK tuning method
%    1   cw_DK          DK            1 (MLE)
%    3   MGK            LT-MGK        [] (not applicable)
%    4   ilrt_MGK       ILRT-MGK      []
%    6   beta_logitMGK  BLT-MGK       []
%    7   gamma_MGK      LGT-MGK       []
%
% Indices 2 and 5 are unused: the numbering is preserved so the result file names
% keep matching those reported in the paper.
% --------------------- Pedro Saa UC 2026 ---------------------------------
specs(1) = mk(1,'cw_DK',1,'DK','population_cwDK','results_bayes_net1');
specs(2) = mk(3,'MGK',[],'LT-MGK','population_MGK','results_bayes_net3');
specs(3) = mk(4,'ilrt_MGK',[],'ILRT-MGK','population_ilrtMGK','results_bayes_net4');
specs(4) = mk(6,'beta_logitMGK',[],'BLT-MGK','population_beta_logitMGK','results_bayes_net6');
specs(5) = mk(7,'gamma_MGK',[],'LGT-MGK','population_gammaMGK','results_bayes_net7');

if (nargin<1)
    spec = specs;
    return
end

hit = find([specs.idx]==kernelIdx,1);
if isempty(hit)
    error('getKernelSpec:unknownIndex',...
          'Unknown kernel index %d. Valid indices are: %s.',kernelIdx,mat2str([specs.idx]));
end
spec = specs(hit);
end

%% ------------------------------------------------------------------------
function s = mk(idx,kernel,method,paperName,varName,resultFile)
s = struct('idx',idx,'kernel',kernel,'method',method,'paperName',paperName,...
           'varName',varName,'resultFile',resultFile);
end

function y = mappingFxn_beta(x)
% Betancourt's forward map from the simplex onto (D-1) independent Beta variates.
%
%   y_1 = 1 - x_1
%   y_j = 1 - x_j / prod_{k<j} y_k        for 1 < j < D
%
% The running product prod_{k<j} y_k is the remaining stick length,
% 1 - sum_{k<j} x_k, which is what makes the y's independent. It must run over
% ALL preceding coordinates: dividing by y_{j-1} alone agrees only up to D = 3.
% --------------------- Pedro Saa UC 2026 ---------------------------------
D      = size(x,2)-1;
y      = zeros(size(x,1),D);
y(:,1) = 1-x(:,1);
for ix = 2:D
    y(:,ix) = 1 - x(:,ix)./prod(y(:,1:ix-1),2);
end
end

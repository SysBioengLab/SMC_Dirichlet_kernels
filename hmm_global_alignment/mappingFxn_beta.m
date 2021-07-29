function y = mappingFxn_beta(x)
D      = size(x,2)-1;
y      = zeros(size(x,1),D);
y(:,1) = 1-x(:,1);
if (D>1)
    for ix = 2:D
        y(:,ix) = 1 - x(:,ix)./prod(y(:,ix-1),2);
    end
end
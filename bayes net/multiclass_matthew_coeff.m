function MCC = multiclass_matthew_coeff(K)
% Compute Matthew's Correlation Coefficient for multiclass classification.
% Inputs:   K confusion matrix (double array)
%
% Outputs:  MCC
%------------------------ Pedro Saa UC 2021 -------------------------------
Ktrans = K';
N  = sum(sum(K));
sum1 = 0; sum2 = 0; sum3 = 0;
for k = 1:size(K,1)
    for l = 1:size(K,2)
        sum1 = sum1+K(k,:)*K(:,l);
        sum2 = sum2+K(k,:)*Ktrans(:,l);
        sum3 = sum3+Ktrans(k,:)*K(:,l);
    end
end
MCC = (N*trace(K)-sum(sum(sum1)))/(sqrt(N^2-sum(sum(sum2)))*sqrt(N^2-sum(sum(sum3))));
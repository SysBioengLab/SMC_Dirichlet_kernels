function p = predictCovidProb(th,off,sym,odo)
% P(COV = positive | symptoms, odours) in closed form, for one parameter vector.
%
% Inputs:   th   1 x 72 flat parameter vector (all Dirichlet blocks)
%           off  node offsets from getCovidOffsets
%           sym  n x 7 symptom states in {1,2}
%           odo  n x 6 odour states in {1,2}
% Output:   p    1 x n predictive probability of a positive RT-qPCR
%
%   P(COV=+|s,o) ~ P(COV=+) prod_k P(s_k|COV=+) sum_m P(sKor=m|COV=+) prod_l P(o_l|sKor=m)
%
% The sum marginalises the unobserved olfactory status: this is what makes the
% likelihood intractable for the sampler, and what makes the prediction invariant
% to the label switching of that ternary node. Verified against BNT's
% junction-tree engine to machine precision (see classifier_performance).
% --------------------- Pedro Saa UC 2026 ---------------------------------
st = off.start; cl = off.cols;
n  = size(sym,1);

pc = [th(st(1)), th(st(1)+1)];                                             % P(COV = 1), P(COV = 2)

A = ones(n,2);                                                             % symptom likelihood per COV state
for k = 1:7
    nd = 1+k;
    for c = 1:2
        p1 = th(st(nd) + (c-1)*cl(nd));                                    % P(symptom = state 1 | COV = c)
        A(:,c) = A(:,c).*(p1*(sym(:,k)==1) + (1-p1)*(sym(:,k)==2));
    end
end

Lo = ones(n,3);                                                            % odour likelihood per sKor state
for l = 1:6
    nd = 9+l;
    for m = 1:3
        p1 = th(st(nd) + (m-1)*cl(nd));                                    % P(odour = state 1 | sKor = m)
        Lo(:,m) = Lo(:,m).*(p1*(odo(:,l)==1) + (1-p1)*(odo(:,l)==2));
    end
end

B = zeros(n,2);                                                            % marginalise the hidden sKor
for c = 1:2
    pk     = th(st(9) + (c-1)*cl(9) + (0:2));
    B(:,c) = Lo*pk(:);
end

num = pc(2)*A(:,2).*B(:,2);
den = pc(1)*A(:,1).*B(:,1) + num;
p   = (num./max(den,realmin))';
end

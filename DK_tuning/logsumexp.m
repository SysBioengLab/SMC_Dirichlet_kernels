function f = logsumexp(x)
% log(sum(exp(x))), shifted by max(x) to avoid overflow.
%--------------------- Pedro Saa UC 2026 ----------------------------------
the_max = max(x);
x       = x - the_max;
f       = the_max + log(sum(exp(x)));
end

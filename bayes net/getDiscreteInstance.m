function discreteInstance = getDiscreteInstance(contInstance)
discreteInstance = zeros(size(contInstance));

% Discretize 1st node
if (contInstance(1)<1+.5)
    discreteInstance(1) = 1;
else
    discreteInstance(1) = 2;
end

% Discretize 2nd node
if (contInstance(2) < 1+2/3)
    discreteInstance(2) = 1;
elseif (contInstance(2) < 1+2*2/3)
    discreteInstance(2) = 2;
else
    discreteInstance(2) = 3;
end

% Discretize 3rd node
if (contInstance(3) < 1+2/3)
    discreteInstance(3) = 1;
elseif (contInstance(3) < 1+2*2/3)
    discreteInstance(3) = 2;
else
    discreteInstance(3) = 3;
end

% Discretize 4th node
if (contInstance(4)<1+.5)
    discreteInstance(4) = 1;
else
    discreteInstance(4) = 2;
end

% Discretize 5th node
if (contInstance(5)<1+.5)
    discreteInstance(5) = 1;
else
    discreteInstance(5) = 2;
end

% Discretize 6th node
if (contInstance(6)<1+.5)
    discreteInstance(6) = 1;
else
    discreteInstance(6) = 2;
end
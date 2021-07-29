function predState = simulateBayesNetCond(probs,condState,N)
% Simulate N instances of the Bayesian net described by the conditional
% probabilities probs
%------------------------- Pedro Saa UC 2021 ------------------------------
unknownStates = find(condState==0);
states        = repmat(condState,N,1);

% Main loop
for ix = 1:N
    
    % If unknown, simulate node 1: ST-stimulant
    if any(unknownStates==1)
        states(ix,1) = sum(rand(1) >= cumsum([0,probs{1}]));
    end
    
    % If unknown, simulate node 2 cond. on node 1: SI-signal
    if any(unknownStates==2)
        if (states(ix,1)==1)
            states(ix,2) = sum(rand(1) >= cumsum([0,probs{2}(1,:)]));       % 1st-row for SI present
            
        elseif (states(ix,1)==2)
            states(ix,2) = sum(rand(1) >= cumsum([0,probs{2}(2,:)]));       % 2nd-row for SI absent
        end
    end
    
    % If unknown, simulate node 3 cond. on node 2: IN-inhibitor
    if any(unknownStates==3)
        if (states(ix,2)==1)
            states(ix,3) = sum(rand(1) >= cumsum([0,probs{3}(1,:)]));       % 1st-row for SI-high
            
        elseif (states(ix,2)==2)
            states(ix,3) = sum(rand(1) >= cumsum([0,probs{3}(2,:)]));       % 2nd-row for SI-medium
            
        elseif (states(ix,2)==3)
            states(ix,3) = sum(rand(1) >= cumsum([0,probs{3}(3,:)]));       % 2nd-row for SI-low
        end
    end
    
    % If unknown, simulate node 4 cond. on node 2&3: RE-receptor
    if any(unknownStates==4)
        if (states(ix,2)==1)&&(states(ix,3)==1)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(1,:)]));       % High-High
            
        elseif (states(ix,2)==1)&&(states(ix,3)==2)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(2,:)]));       % High-Med
            
        elseif (states(ix,2)==1)&&(states(ix,3)==3)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(3,:)]));       % High-Low
            
        elseif (states(ix,2)==2)&&(states(ix,3)==1)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(4,:)]));       % Med-High
            
        elseif (states(ix,2)==2)&&(states(ix,3)==2)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(5,:)]));       % Med-Med
            
        elseif (states(ix,2)==2)&&(states(ix,3)==3)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(6,:)]));       % Med-Low
            
        elseif (states(ix,2)==3)&&(states(ix,3)==1)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(7,:)]));       % Low-High
            
        elseif (states(ix,2)==3)&&(states(ix,3)==2)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(8,:)]));       % Low-Med
            
        elseif (states(ix,2)==3)&&(states(ix,3)==3)
            states(ix,4) = sum(rand(1) >= cumsum([0,probs{4}(9,:)]));       % Low-Low
        end
    end
    
    % If unknown, simulate node 5 cond. on node 4: GP-G protein
    if any(unknownStates==5)
        if (states(ix,4)==1)
            states(ix,5) = sum(rand(1) >= cumsum([0,probs{5}(1,:)]));       % 1st-row for GP active
            
        elseif (states(ix,4)==2)
            states(ix,5) = sum(rand(1) >= cumsum([0,probs{5}(2,:)]));       % 2nd-row for GP not active
        end
    end
    
    % If unknown, simulate node 6 cond. on node 5: CR-cellular response
    if any(unknownStates==6)
        if (states(ix,5)==1)
            states(ix,6) = sum(rand(1) >= cumsum([0,probs{6}(1,:)]));       % 1st-row for GP-P bound (YES)
            
        elseif (states(ix,5)==2)
            states(ix,6) = sum(rand(1) >= cumsum([0,probs{6}(2,:)]));       % 2nd-row for GP-P not bound (NO)
        end
    end
end

% Return predicted state (expected value)
predState = mean(states);

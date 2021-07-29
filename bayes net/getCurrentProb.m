function logProb = getCurrentProb(prevStates,newState,logProb,probs,level)
% Simulate N instances of the Bayesian net described by the conditional
% prbabilities probs
%------------------------- Pedro Saa UC 2021 ------------------------------

% Node 1: ST-stimulant
if (level==1)    
    logProb  = log(probs{1}(newState)) + logProb;
   
% Node 2 cond. on node 1: SI-signal
elseif (level==2)
    if (prevStates(1)==1)
        logProb = log(probs{2}(1,newState)) + logProb;
                
    elseif (prevStates(1)==2)        
        logProb = log(probs{2}(2,newState)) + logProb;
    end
    
% Node 3 cond. on node 2: IN-inhibitor
elseif (level==3)
    if (prevStates(2)==1)
        logProb = log(probs{3}(1,newState)) + logProb;
        
    elseif (prevStates(2)==2)
        logProb = log(probs{3}(2,newState)) + logProb;
        
    elseif (prevStates(2)==3)
        logProb = log(probs{3}(3,newState)) + logProb;
    end
    
% Node 4 cond. on node 2&3: RE-receptor
elseif (level==4)   
    if (prevStates(2)==1)&&(prevStates(3)==1)        
        logProb = log(probs{4}(1,newState)) + logProb;
        
    elseif (prevStates(2)==1)&&(prevStates(3)==2)
        logProb = log(probs{4}(2,newState)) + logProb;
        
    elseif (prevStates(2)==1)&&(prevStates(3)==3)
        logProb = log(probs{4}(3,newState)) + logProb;
        
    elseif (prevStates(2)==2)&&(prevStates(3)==1)
        logProb = log(probs{4}(4,newState)) + logProb;
        
    elseif (prevStates(2)==2)&&(prevStates(3)==2)
        logProb = log(probs{4}(5,newState)) + logProb;
        
    elseif (prevStates(2)==2)&&(prevStates(3)==3)
        logProb = log(probs{4}(6,newState)) + logProb;
        
    elseif (prevStates(2)==3)&&(prevStates(3)==1)
        logProb = log(probs{4}(7,newState)) + logProb;
        
    elseif (prevStates(2)==3)&&(prevStates(3)==2)
        logProb = log(probs{4}(8,newState)) + logProb;
        
    elseif (prevStates(2)==3)&&(prevStates(3)==3)
        logProb = log(probs{4}(9,newState)) + logProb;
    end

% Node 5 cond. on node 4: GP-G protein
elseif (level==5)    
    if (prevStates(4)==1)
        logProb = log(probs{5}(1,newState)) + logProb;
        
    elseif (prevStates(4)==2)
        logProb = log(probs{5}(2,newState)) + logProb;
    end

% Node 6 cond. on node 5: CR-cellular response
elseif (level==6)
    if (prevStates(5)==1)
        logProb = log(probs{6}(1,newState)) + logProb;
        
    elseif (prevStates(5)==2)
        logProb = log(probs{6}(2,newState)) + logProb;
    end
end
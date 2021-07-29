function [probState,maxProb] = posteriorProbBayesNet(probs,state)
% Computes the maximum posterior probability given a potential state
% Input:    prob parameter (cell array)
% Output:   post prob
%--------------------- Pedro Saa UC 2021 ----------------------------------
% Definition of global variables
postProbs = [];
allStates = [];

% Search parameters
depth     = numel(state);
currState = zeros(1,depth);

% Define node space
nodeSpace{1} = [1,2];
nodeSpace{2} = [1,2,3];
nodeSpace{3} = [1,2,3];
nodeSpace{4} = [1,2];
nodeSpace{5} = [1,2];
nodeSpace{6} = [1,2];

% Perform depth-first search: recursive call
depthSearch(currState,0,0);

% Return maximum posterior prob
[maxProb,ix] = max(postProbs);
probState    = allStates(ix,:);

% Depth-first search
    function depthSearch(currState,logProb,level)

        % Solve LP problem at the bottom node
        if (depth == level)

            % Check feasibility and save if tfm feasible
            allStates = [allStates;currState];
            postProbs = [postProbs;logProb];
            
            % Recursive call
        else
            % Increase a level
            level = level+1;
            
            % Check if the state is unknown (if so, branch out)
            if (state(level)==0)                
                for jx = 1:numel(nodeSpace{level})
                    
                    %Assign new state
                    newState = nodeSpace{level}(jx);
                    
                    % Compute cumulative probability
                    logProbTemp = getCurrentProb(currState,newState,logProb,probs,level);
                    
                    % Update current state
                    currState(level) = newState;
                    
                    % Advance a level in the recursion
                    depthSearch(currState,logProbTemp,level);
                end
            
            % State is known
            else
                
                %Assign new state
                newState = nodeSpace{level}(state(level));
                
                % Compute cumulative probability
                logProb  = getCurrentProb(currState,newState,logProb,probs,level);
                
                % Update current state
                currState(level) = newState;
                
                % Advance a level in the recursion
                depthSearch(currState,logProb,level);
            end
        end
    end
end
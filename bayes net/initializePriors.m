function priors = initializePriors(fullData)
% Initialize priors for Bayesian network based on complete data
%------------------------- Pedro Saa UC 2021 ------------------------------
priors{1} = [sum(fullData(:,1)==1)+1,sum(fullData(:,1)==2)+1];             % STIMULANT

priors{2} = [sum((fullData(:,1)==1)&(fullData(:,2)==1))+1,...              % SIGNAL (ST present)
             sum((fullData(:,1)==1)&(fullData(:,2)==2))+1,...
             sum((fullData(:,1)==1)&(fullData(:,2)==3))+1];                
priors{2} = [priors{2};                                                    % SIGNAL (ST absent)
             sum((fullData(:,1)==2)&(fullData(:,2)==1))+1,...              
             sum((fullData(:,1)==2)&(fullData(:,2)==2))+1,...
             sum((fullData(:,1)==2)&(fullData(:,2)==3))+1];                

priors{3} = [sum((fullData(:,2)==1)&(fullData(:,3)==1))+1,...              % INHIBITOR (SIGNAL high)
             sum((fullData(:,2)==1)&(fullData(:,3)==2))+1,...
             sum((fullData(:,2)==1)&(fullData(:,3)==3))+1];                
priors{3} = [priors{3};...                                                 % INHIBITOR (SIGNAL medium)
             sum((fullData(:,2)==2)&(fullData(:,3)==1))+1,...
             sum((fullData(:,2)==2)&(fullData(:,3)==2))+1,...
             sum((fullData(:,2)==2)&(fullData(:,3)==3))+1];                
priors{3} = [priors{3};...                                                 % INHIBITOR (SIGNAL low)
             sum((fullData(:,2)==3)&(fullData(:,3)==1))+1,...
             sum((fullData(:,2)==3)&(fullData(:,3)==2))+1,...
             sum((fullData(:,2)==3)&(fullData(:,3)==3))+1];

priors{4} = [sum((fullData(:,2)==1)&(fullData(:,3)==1)&(fullData(:,4)==1))+1,...    % RECEPTOR (SIGNAL high INHIBITOR high)
             sum((fullData(:,2)==1)&(fullData(:,3)==1)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...
             sum((fullData(:,2)==1)&(fullData(:,3)==2)&(fullData(:,4)==1))+1,...    % RECEPTOR (SIGNAL high INHIBITOR medium)
             sum((fullData(:,2)==1)&(fullData(:,3)==2)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...
             sum((fullData(:,2)==1)&(fullData(:,3)==3)&(fullData(:,4)==1))+1,...    % RECEPTOR (SIGNAL high INHIBITOR low)
             sum((fullData(:,2)==1)&(fullData(:,3)==3)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL med INHIBITOR high)
             sum((fullData(:,2)==2)&(fullData(:,3)==1)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==2)&(fullData(:,3)==1)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL med INHIBITOR med)
             sum((fullData(:,2)==2)&(fullData(:,3)==2)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==2)&(fullData(:,3)==2)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL med INHIBITOR low)
             sum((fullData(:,2)==2)&(fullData(:,3)==3)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==2)&(fullData(:,3)==3)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL low INHIBITOR high)
             sum((fullData(:,2)==3)&(fullData(:,3)==1)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==3)&(fullData(:,3)==1)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL low INHIBITOR med)
             sum((fullData(:,2)==3)&(fullData(:,3)==2)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==3)&(fullData(:,3)==2)&(fullData(:,4)==2))+1];
priors{4} = [priors{4};...                                                          % RECEPTOR (SIGNAL low INHIBITOR low)
             sum((fullData(:,2)==3)&(fullData(:,3)==3)&(fullData(:,4)==1))+1,...
             sum((fullData(:,2)==3)&(fullData(:,3)==3)&(fullData(:,4)==2))+1];

priors{5} = [sum((fullData(:,4)==1)&(fullData(:,5)==1))+1,...                       % G PROTEIN (RECEPTIOR active)
             sum((fullData(:,4)==1)&(fullData(:,5)==2))+1];
priors{5} = [priors{5};...                                                          % G PROTEIN (RECEPTIOR inactive)
             sum((fullData(:,4)==2)&(fullData(:,5)==1))+1,...                       
             sum((fullData(:,4)==2)&(fullData(:,5)==2))+1];
         
priors{6} = [sum((fullData(:,5)==1)&(fullData(:,6)==1))+1,...                       % CELLULAR RESPONSE (G PROTEIN binds)
             sum((fullData(:,5)==1)&(fullData(:,6)==2))+1];
priors{6} = [priors{6};...                                                          % CELLULAR RESPONSE (G PROTEIN doesn't bind)
             sum((fullData(:,5)==2)&(fullData(:,6)==1))+1,...                       
             sum((fullData(:,5)==2)&(fullData(:,6)==2))+1];
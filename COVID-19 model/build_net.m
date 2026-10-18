% Build the COVID-19 screening Bayes net STRUCTURE and save it as Bayes_net.mat.
% The CPD tables are assigned later by run_bayes_net.m and bayes_net_example.m.
%------------------------- Pedro Saa UC 2026 ------------------------------

% Add the relevant folder to the path
addpath('./bnt-master/')                                                     % exposes genpathKPM
addpath(genpathKPM('./bnt-master/'))                                         % exposes BNT and its subfolders (skips Old/ and CVS/)

N = 15;
dag = zeros(N,N);

% Name nodes
COV  = 1;
sCol = 2;
sCou = 3;
sFev = 4;
sMus = 5;
sBre = 6;
sAg  = 7;
sAn  = 8;
sKor = 9;
oB   = 10;
oC   = 11;
oM   = 12;
oO   = 13;
oP   = 14;
oV   = 15;

% Adjacency matrix
dag(COV,[sCol,sCou,sFev,sMus,sBre,sAg,sAn,sKor]) = 1;
dag(sKor,[oB,oC,oM,oO,oP,oV]) = 1;

% Type of nodes. All dichotomous except the hidden sKor (objective olfactory
% status): anosmic / hyposmic / normosmic.
discrete_nodes = 1:N;
node_sizes     = 2*ones(1,N);
node_sizes(sKor) = 3;

% Build Bayes net
bnet = mk_bnet(dag, node_sizes, 'names', {'COV','sCol','sCou',...
    'sFev','sMus','sBre','sAg','sAn','sKor','oB','oC','oM','oO','oP','oV'},...
    'discrete', discrete_nodes);

% To visualise the network structure: draw_graph(bnet.dag)
save Bayes_net
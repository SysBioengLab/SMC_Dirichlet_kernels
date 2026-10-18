function make_network_figure(outFile)
% Draw the structure of the COVID-19 screening Bayes net.
%
% Input:  outFile  destination .png (default ./results/network_structure.png)
%
% A labelled reference drawing; it is not itself a figure of the paper. The
% network appears in the paper as panel A of the Bayes-net figure, drawn into
% the final artwork by hand.
%
% The DAG is the one built by build_net.m: COV is the root of the seven
% self-reported symptoms and of the olfactory status sKor, itself the parent of
% the six odour outcomes. Every node is dichotomous except sKor, which has three
% states and is never observed.
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1) || isempty(outFile); outFile = './results/network_structure.png'; end

names = {'COV','sCol','sCou','sFev','sMus','sBre','sAg','sAn','sKor',...
         'oB','oC','oM','oO','oP','oV'};
sym   = 2:8;                                                               % self-reported symptoms
odo   = 10:15;                                                             % psychophysical odour outcomes
kor   = 9;

s = [ones(1,numel(sym)), 1,            kor*ones(1,numel(odo))];
t = [sym,                kor,          odo];
G = digraph(s,t,[],names);

% Layered layout, laid out by hand so the two fans do not overlap
x = zeros(1,15); y = zeros(1,15);
x(1)   = 3.0;              y(1)   = 2;                                     % COV
x(sym) = linspace(0,6,7);  y(sym) = 1;                                      % symptoms
x(kor) = 8.0;              y(kor) = 1;                                      % hidden node
x(odo) = linspace(5.5,10.5,6); y(odo) = 0;                                  % odours

f = figure('Name','COVID-19 Bayes net structure','Position',[80 80 1100 460],'Color','w');
p = plot(G,'XData',x,'YData',y,'ArrowSize',9,'LineWidth',1.1,...
         'EdgeColor',[.45 .45 .45],'NodeFontSize',10,'MarkerSize',7);
highlight(p,1,'NodeColor',[.85 .33 .27]);                                  % observed root
highlight(p,sym,'NodeColor',[.35 .45 .70]);
highlight(p,odo,'NodeColor',[.35 .45 .70]);
highlight(p,kor,'NodeColor',[.95 .70 .20],'MarkerSize',11);                % hidden, ternary
axis off; axis([-0.6 11.2 -0.55 2.45])

text(3.0,2.32,'observed status (RT-qPCR)','HorizontalAlignment','center','FontSize',8,'FontAngle','italic');
text(8.0,1.32,'hidden, 3 states','HorizontalAlignment','center','FontSize',8,'FontAngle','italic');
text(3.0,-0.42,'self-reported symptoms','HorizontalAlignment','center','FontSize',8,'FontAngle','italic');
text(8.0,-0.42,'psychophysical odour outcomes','HorizontalAlignment','center','FontSize',8,'FontAngle','italic');

d = fileparts(outFile);
if ~isempty(d) && ~isfolder(d); mkdir(d); end
print(f,outFile,'-dpng','-r200');
savefig(f,strrep(outFile,'.png','.fig'));
fprintf('Wrote %s\n',outFile);
end

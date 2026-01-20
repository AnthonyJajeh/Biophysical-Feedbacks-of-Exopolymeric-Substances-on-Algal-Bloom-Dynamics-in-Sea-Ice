clear all;clc;close all;

fig1 = openfig('SIa.fig','visible');
fig2 = openfig('SIb.fig','visible');
fig3 = openfig('SIc.fig','visible');

% Make it current
figure(fig3);

% Get axes
ax = findall(fig3, 'Type', 'Axes');
ax = ax(1);

% DELETE any previous star(s) we added
delete(findall(ax, 'Type', 'Line', 'Tag', 'HopfStar'));


nice_graphing('SIa',fig1);
nice_graphing('SIb',fig2);
nice_graphing('SIc',fig3);



function nice_graphing(fname, fig)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .8; % feel free to play with this ratio
set(findall(fig,'-property','FontSize'),'FontSize',24) % adjust fontsize to your document
set(findall(fig,'-property','Box'),'Box','on') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
lgd = findall(fig, 'Type', 'Legend');
set(lgd, 'Box', 'off');     % ensure no border if a legend exists
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(fig, strcat(fname,'.png'), 'Resolution', 300);
exportgraphics(fig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(fig, strcat(fname,'.fig'));
end

% exportgraphics(fig1, strcat('SIa','.png'), 'Resolution', 300);
% exportgraphics(fig1, strcat('SIa','.pdf'), 'ContentType', 'vector');
% saveas(fig1, strcat('SIa','.fig'));
% 
% 
% exportgraphics(fig2, strcat('SIb','.png'), 'Resolution', 300);
% exportgraphics(fig2, strcat('SIb','.pdf'), 'ContentType', 'vector');
% saveas(fig2, strcat('SIb','.fig'));
% 
% 
% exportgraphics(fig3, strcat('SIc','.png'), 'Resolution', 300);
% exportgraphics(fig3, strcat('SIc','.pdf'), 'ContentType', 'vector');
% saveas(fig3, strcat('SIc','.fig'));
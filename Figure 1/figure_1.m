%Author:Anthony Jajeh 
%Date: Aug 18, 2025
%Figure 1 by using MATLAB built
%in ODE solver
clear all;clc;close all;



% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 

%Parameter values for fig 1b
% phi = .01;
% psi = .01;
% mu = .001;
% gamma = .01; 
% nu = .2; 
% rho = .75; 
% xi = .2;
% delta = .007; 
% eta = .03;

% 
% % %Parameter values for fig 1a
phi = .0001;
psi = .05;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;% 



%nondimensional conversion values 
a = phi/(gamma*delta);
b = psi/delta;
c = nu/delta;
d = (rho*gamma)/(mu*eta);
f = xi * c;
epsilon = eta/delta;

%domain 
n=10;
domain = [0 n];

%Initial conditions
IC_N = .2/gamma;
IC_A = .03/gamma;
IC_E = .8;

%initial condition vector 
IC_fast = [IC_N IC_A];
x = (f*(a*f-a-b))/(c*(f-1));
%scaled values
%Solving QSS fast-model 
[IVsol_fast, DVsol_fast] = ode23s(@(t, y) DEdef_fast(t, y, a,b,c,f,IC_E), domain, IC_fast);
N_sol_fast = DVsol_fast(:, 1);
A_sol_fast = DVsol_fast(:, 2);
IVsol_fast = IVsol_fast;


% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_fast, N_sol_fast, 'color', nutrientcolordet, 'linewidth', 4);
ylim([0, max(N_sol_fast) * 1.2]);
ylabel('nutrients','Color','k');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black

% Plot algae on the right y-axis
yyaxis right;
plot(IVsol_fast, A_sol_fast, 'color', algaecolordet, 'linewidth', 4);
ylabel('algae','Color','k');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black

% Set common properties
xlim([0, n]);
xlabel('time');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Add legend
legend('Nutrients', 'Algae', 'Location', 'northeast','Box','off');


hold off;
% Set the figure size and save as PDF, PNG, and FIG
set(gcf, 'units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]); % Set the paper to have width 6 and height 4
fname = 'fig1b';
nice_graphing(fname, fig)

function nice_graphing(fname, fig)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .8; % feel free to play with this ratio
set(findall(fig,'-property','FontSize'),'FontSize',24) % adjust fontsize to your document
set(findall(fig,'-property','Box'),'Box','on') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
 lgd = findall(fig, 'Type', 'Legend');
    set(lgd, 'Box', 'off');     % ensure no border if a legend exists
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(fig, strcat(fname,'.png'), 'ContentType', 'vector');
end

%Defining QSS-NA-model
function [Dode] = DEdef_fast(I,D,a,b,c,f,E_0)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);

%set of odes
dNdt = a*exp(-E_0)-(c*N*A)/(N+1) - b*N*exp(-E_0);
dAdt = (f*N*A)/(N+1)-A;

% odes in vector form
Dode = [dNdt; dAdt];
end

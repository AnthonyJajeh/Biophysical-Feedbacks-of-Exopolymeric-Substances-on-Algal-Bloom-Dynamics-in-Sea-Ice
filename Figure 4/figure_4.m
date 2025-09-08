%Author:Anthony Jajeh
%Date: Febuary 5th, 2025
% Steady state analysis of NAE-model using exponential function as
% inflow/outflow of nutrients functions
clear all; clc; close all;
n=300;

algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS
% 
% %Parameter values for fig 4a
phi = .0001;
psi = .1;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;
% % % 

% %Parameter values for fig 4b
% phi = .0001;
% psi = .01;
% mu = .001;
% gamma = .01; 
% nu = .2; 
% rho = .75; 
% xi = .2;
% delta = .007; 
% eta = .03;

% %Parameter values for fig 4c
% phi = .01;
% psi = .01;
% mu = .001;
% gamma = .01; 
% nu = .2; 
% rho = .75; 
% xi = .2;
% delta = .007; 
% eta = .03;

%nondimensional conversion values 
a = phi/(gamma*delta);
b = psi/delta;
c = nu/delta;
d = (rho*gamma)/(mu*eta);
f = xi * c;
epsilon = eta/delta;

domain = [0 n];
%Initial conditions
IC_N = .2/gamma;
IC_A = .03/gamma;
IC_E = .8/mu;

IC_exp = [IC_N IC_A IC_E];
%calculating NAE-model solution plots 
[IVsol_exp, DVsol_exp] = ode23s(@(t, y) DEdef_exp(t, y, a,b,c,f,d,epsilon), domain, IC_exp);
N_sol_exp = DVsol_exp(:, 1);
A_sol_exp = DVsol_exp(:, 2);
E_sol_exp = DVsol_exp(:, 3);
IVsol_exp = IVsol_exp;


%plotting solution curves of full-model 
% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp, N_sol_exp, 'color', nutrientcolordet, 'linewidth', 3);
plot(IVsol_exp, A_sol_exp, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([N_sol_exp;A_sol_exp]) * 1.2]);
ylabel('nutrients \& algae','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp, E_sol_exp, 'color', EPScolordet, 'linewidth', 3);
ylim([0, max(E_sol_exp) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('EPS','FontSize',20,'Color','k');

% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors

% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)



%Defining NAE-model



fname = 'fig4a';
nice_graphing(fname,fig)


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
exportgraphics(fig, strcat(fname,'.png'), 'ContentType', 'vector');
end

function [Dode] = DEdef_exp(I,D,a,b,c,f,d,epsilon)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = (1/epsilon)* (a*exp(-E)-(c*A*N)/(N+1)-b*N*exp(-E));
dAdt = (1/epsilon)* ((f*N*A)/(N+1)-A);
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
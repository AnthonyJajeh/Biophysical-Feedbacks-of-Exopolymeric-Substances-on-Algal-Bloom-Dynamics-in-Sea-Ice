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
%Parameter values for fig 3a
phi_0 = .0001;
psi_0 = .05;
mu_0 = .001;
gamma_0 = .01; 
nu_0 = .2; 
rho_0 = .75; 
xi_0 = .2;
delta_0 = .007; 
eta_0 = .03;
% % 

%Parameter values for fig 3b
phi_1 = .0001;
psi_1 = .001;
mu_1 = .001;
gamma_1 = .01; 
nu_1 = .2; 
rho_1 = .75; 
xi_1 = .2;
delta_1 = .007; 
eta_1 = .03;

% Parameter values for fig 3c
phi_2 = .01;
psi_2 = .01;
mu_2 = .001;
gamma_2 = .01; 
nu_2 = .2; 
rho_2 = .75; 
xi_2 = .2;
delta_2 = .007; 
eta_2 = .03;

%nondimensional conversion values 
a_0 = phi_0/(gamma_0*delta_0);
b_0 = psi_0/delta_0;
c_0 = nu_0/delta_0;
d_0 = (rho_0*gamma_0)/(mu_0*eta_0);
f_0 = xi_0 * c_0;
epsilon_0 = eta_0/delta_0;

%nondimensional conversion values 
a_1 = phi_1/(gamma_1*delta_1);
b_1 = psi_1/delta_1;
c_1 = nu_1/delta_1;
d_1 = (rho_1*gamma_1)/(mu_1*eta_1);
f_1 = xi_1 * c_1;
epsilon_1 = eta_1/delta_1;

%nondimensional conversion values 
a_2 = phi_2/(gamma_2*delta_2);
b_2 = psi_2/delta_2;
c_2 = nu_2/delta_2;
d_2 = (rho_2*gamma_2)/(mu_2*eta_2);
f_2 = xi_2 * c_2;
epsilon_2 = eta_2/delta_2;

domain = [0 n];
%Initial conditions
IC_N = .2/gamma_0;
IC_A = .0002/gamma_0;
IC_E = .002/mu_0;

IC_exp = [IC_N IC_A IC_E];
%calculating NAE-model solution plots a
[IVsol_exp_triv, DVsol_exp_triv] = ode23s(@(t, y) DEdef_exp(t, y, a_0,b_0,c_0,f_0,d_0,epsilon_0), domain, IC_exp);
N_sol_exp_triv = DVsol_exp_triv(:, 1);
A_sol_exp_triv = DVsol_exp_triv(:, 2);
E_sol_exp_triv = DVsol_exp_triv(:, 3);

%calculating NAE-model solution plots b
[IVsol_exp_sp, DVsol_exp_sp] = ode23s(@(t, y) DEdef_exp(t, y, a_1,b_1,c_1,f_1,d_1,epsilon_1), domain, IC_exp);
N_sol_exp_sp = DVsol_exp_sp(:, 1);
A_sol_exp_sp = DVsol_exp_sp(:, 2);
E_sol_exp_sp = DVsol_exp_sp(:, 3);

%calculating NAE-model solution plots c
[IVsol_exp_hopf, DVsol_exp_hopf] = ode23s(@(t, y) DEdef_exp(t, y, a_2,b_2,c_2,f_2,d_2,epsilon_2), domain, IC_exp);
N_sol_exp_hopf = DVsol_exp_hopf(:, 1);
A_sol_exp_hopf = DVsol_exp_hopf(:, 2);
E_sol_exp_hopf = DVsol_exp_hopf(:, 3);



%plotting solution curves of full-model 
% Create a new figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_triv = figure;
set(fig_triv, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp_triv, N_sol_exp_triv, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, .5]);
ylabel('nutrients','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp_triv, A_sol_exp_triv, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_exp_triv, E_sol_exp_triv, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0,1]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');


% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
% %Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_sp = figure;
set(fig_sp, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp_sp, N_sol_exp_sp, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, 1]);
ylabel('nutrients','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp_sp, A_sol_exp_sp, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_exp_sp, E_sol_exp_sp, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0,6]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');


% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
% %Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_hopf = figure;
set(fig_hopf, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp_hopf, N_sol_exp_hopf, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, 1]);
ylabel('nutrients','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp_hopf, A_sol_exp_hopf, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_exp_hopf, E_sol_exp_hopf, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0,18]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');


% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
% %Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors


fname1 = 'fig3a';
fname2 = 'fig3b';
fname3 = 'fig3c';
nice_graphing(fname1,fig_triv)
nice_graphing(fname2,fig_sp)
nice_graphing(fname3,fig_hopf)


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
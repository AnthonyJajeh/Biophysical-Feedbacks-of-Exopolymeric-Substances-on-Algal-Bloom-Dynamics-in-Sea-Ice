%Author: Anthony Jajeh
%Date: June 3rd, 2025
% Creates figure 9
% Steady state analysis of realistic-model using exponential function as
% inflow/outflow of nutrients functions
clear all; clc; close all;

n_0 = 365;
n_1 = 3000;
domain_0 = [0 n_0];
domain_1 = [0 n_1];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS

%Parameter values for fig 6a
phi = .01;
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;


%Initial conditions
IC_N = .2;
IC_A = .0002;
IC_E = .002;


IC_exp = [IC_N IC_A IC_E];

%calculating realistic-model solution plots 
[IVsol_exp_short, DVsol_exp_short] = ode23s(@(t, y) DEdef_realistic(t, y, phi, psi ,mu, gamma, nu, rho, xi, delta, eta), domain_0, IC_exp);
N_sol_exp_short = DVsol_exp_short(:, 1);
A_sol_exp_short = DVsol_exp_short(:, 2);
E_sol_exp_short = DVsol_exp_short(:, 3);

[IVsol_exp_long, DVsol_exp_long] = ode23s(@(t, y) DEdef_realistic(t, y, phi, psi ,mu, gamma, nu, rho, xi, delta, eta), domain_1, IC_exp);
N_sol_exp_long = DVsol_exp_long(:, 1);
A_sol_exp_long = DVsol_exp_long(:, 2);
E_sol_exp_long = DVsol_exp_long(:, 3);

%plotting solution curves of realistic-model 
% Create a new figure
fig_short = figure;
set(fig_short, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp_short, N_sol_exp_short, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_exp_short)*1.2]);
ylabel('nutrients (mg N/L)','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
hold on;
plot(IVsol_exp_short, A_sol_exp_short, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_exp_short, E_sol_exp_short, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max(max(A_sol_exp_short),max(E_sol_exp_short))*1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae (mg chl a/L)\ \&\ EPS (mg XG/L)}$','FontSize',20,'Color','k','Interpreter','latex');

% Set common properties
xlim([0, n_0]);
xlabel('time','FontSize',20,'Color','k');
%Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors

% Create a new figure
fig_long= figure;
set(fig_long, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp_long, N_sol_exp_long, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_exp_short)*1.2]);
ylabel('nutrients (mg N/L)','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
hold on;
plot(IVsol_exp_long, A_sol_exp_long, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_exp_long, E_sol_exp_long, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max(max(A_sol_exp_long),max(E_sol_exp_long))*1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae (mg chl a/L)\ \&\ EPS (mg XG/L)}$','FontSize',20,'Color','k','Interpreter','latex');

% Set common properties
xlim([0, n_1]);
xlabel('time','FontSize',20,'Color','k');
%Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colorsfname = 'fig8b';

fname1 = 'fig6a';
fname2 = 'fig6b';
nice_graphing(fname1, fig_short)
nice_graphing(fname2, fig_long)


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
exportgraphics(fig, strcat(fname,'.png'), 'Resolution', 300);
exportgraphics(fig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(fig, strcat(fname,'.fig'));

end
%Defining realistic-full-model
function [Dode] = DEdef_realistic(I,D,phi,psi,mu, gamma, nu, rho, xi, delta, eta)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = phi * exp(-E/mu) - (nu*N*A)/(N+gamma) - psi*N*exp(-E/mu);
dAdt = (xi*nu*A*N)/(N+gamma) - delta*A;
dEdt = rho*A - eta*E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
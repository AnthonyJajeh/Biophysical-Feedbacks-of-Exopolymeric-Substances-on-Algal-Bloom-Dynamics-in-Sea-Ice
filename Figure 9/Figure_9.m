%Author: Anthony Jajeh
%Date: June 3rd, 2025
% Creates figure 9
% Steady state analysis of realistic-model using exponential function as
% inflow/outflow of nutrients functions
clear all; clc; close all;
n=250;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS

%Parameter values
phi = .05;
psi = .5;
mu = .4;
gamma = .01; 
nu_1 = .24; 
nu_2 = .026; 
xi = .15;
delta = .007; 
eta = .03;


%creating initial condition vectors
IC_N = .165;
IC_A = .0225;
IC_E = .79;
IC_exp = [IC_N IC_A IC_E];

%calculating realistic-model solution plots 
[IVsol_exp, DVsol_exp] = ode23(@(t, y) DEdef_realistic(t, y, phi, psi ,mu, gamma, nu_1, nu_2, xi, delta, eta), domain, IC_exp);
N_sol_exp = DVsol_exp(:, 1);
A_sol_exp = DVsol_exp(:, 2);
E_sol_exp = DVsol_exp(:, 3);

%plotting solution curves of realistic-model 
% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp, N_sol_exp, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_exp) * 1.2]);
ylabel('nutrients','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp, A_sol_exp, 'color', algaecolordet, 'linewidth', 3);
hold on;
plot(IVsol_exp, E_sol_exp, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([A_sol_exp; E_sol_exp]) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('algae & EPS','FontSize',20,'Color','k');

% Set common properties
xlim([0, n]);
xlabel('time (days)','FontSize',20,'Color','k');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors


% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

%Defining realistic-full-model
function [Dode] = DEdef_realistic(I,D,phi,psi,mu, gamma, nu_1, nu_2, xi, delta, eta)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = phi * exp(-E/mu) - (nu_1*N*A)/(N+gamma) - psi*N*exp(-E/mu);
dAdt = (xi*nu_1*A*N)/(N+gamma) - delta*A;
dEdt = nu_2*A - eta*E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
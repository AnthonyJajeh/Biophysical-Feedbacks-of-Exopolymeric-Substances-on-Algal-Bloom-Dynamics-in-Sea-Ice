%Author:Anthony Jajeh 
%Date: Febuary 5th, 2024
%Quasi-steady-state figures by using MATLAB built
%in ODE solver
clear all;clc;close all;

%domain 
n=25;
domain = [0 n];

% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 

%parameters values 
a = 2;
b = .2;
c = 1.3;
c_p = 2;
%Initial conditions
IC_N = 5;
IC_A = 1;
IC_E = 0;

%initial condition vector 
IC_fast = [IC_N IC_A];

%Solving QSS NA-model 
[IVsol_fast, DVsol_fast] = ode23(@(t, y) DEdef_fast(t, y, a,b,c,c_p,IC_E), domain, IC_fast);
N_sol_fast = DVsol_fast(:, 1);
A_sol_fast = DVsol_fast(:, 2);

% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_fast, N_sol_fast, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_fast) * 1.2]);
ylabel('nutrients','Color','k');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black

% Plot algae on the right y-axis
yyaxis right;
plot(IVsol_fast, A_sol_fast, 'color', algaecolordet, 'linewidth', 3);
ylabel('algae','Color','k');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black

% Set common properties
xlim([0, n]);
xlabel('time (days)');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Add legend
legend('Nutrients', 'Algae', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)


hold off;
% Set the figure size and save as PDF, PNG, and FIG
set(gcf, 'units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]); % Set the paper to have width 6 and height 4

%Defining QSS-NA-model
function [Dode] = DEdef_fast(I,D,a,b,c,c_p,E_0)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);

%set of odes
dNdt = a*exp(E_0)-(c*N*A)/(N+1) - b*N*exp(E_0);
dAdt = (c_p*N*A)/(N+1)-A;

% odes in vector form
Dode = [dNdt; dAdt];
end

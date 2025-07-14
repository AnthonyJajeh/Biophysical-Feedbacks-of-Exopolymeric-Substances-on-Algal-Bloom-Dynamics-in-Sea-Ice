%Author:Anthony Jajeh 
%Date: June 3rd, 2025
%Figure 1 by using MATLAB built
%in ODE solver
clear all;clc;close all;

%domain 
n=15;
domain = [0 n];

% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 

%Parameter values for fig 1a
phi = .0001;
psi = .1;
mu = .08;
gamma = .01; 
nu_1 = .2; 
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;

% % 
% % %Parameter values for fig 1b
% phi = .0001;
% psi = .001;
% mu = .15;
% gamma = .01; 
% nu_1 = .2; 
% nu_2 = .05; 
% xi = .2;
% delta = .007; 
% eta = .03;% 

%nondimensional conversion values 
a = phi/(gamma*delta);
b = psi/delta;
c = nu_1/delta;
d = (nu_2*gamma)/(mu*eta);
f = xi * c;

%Initial conditions
IC_N = 5;
IC_A = .03;
IC_E = .01;

value = IC_A*((f*IC_N)/(IC_N+1)-1);
fprintf('N = %.2f mg N/L\n', value);
%initial condition vector 
IC_fast = [IC_N IC_A];

%scaled values
%Solving QSS fast-model 
[IVsol_fast, DVsol_fast] = ode45(@(t, y) DEdef_fast(t, y, a,b,c,f,IC_E), domain, IC_fast);
N_sol_fast = DVsol_fast(:, 1)*gamma;
A_sol_fast = DVsol_fast(:, 2)*gamma;

% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_fast, N_sol_fast, 'color', nutrientcolordet, 'linewidth', 4);
ylim([0, max(N_sol_fast) * 1.2]);
ylabel('nutrients (mg N/L)','Color','k');
set(gca, 'fontsize', 18, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black

% Plot algae on the right y-axis
yyaxis right;
plot(IVsol_fast, A_sol_fast, 'color', algaecolordet, 'linewidth', 4);
ylabel('algae (mg chl-a/L)','Color','k');
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
function [Dode] = DEdef_fast(I,D,a,b,c,f,E_0)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);

%set of odes
dNdt = a*exp(E_0)-(c*N*A)/(N+1) - b*N*exp(E_0);
dAdt = (f*N*A)/(N+1)-A;

% odes in vector form
Dode = [dNdt; dAdt];
end

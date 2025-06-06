%Author:Anthony Jajeh
%Date: June 5th, 2025
% Steady state analysis of full-model using hill function as
% inflow/outflow of nutrients functions
clear all; clc; close all;
n=2500;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 

%Parameter values 
phi = .0015;
psi = .001;
mu = .000088;
gamma = .01; 
nu_1 = .2; 
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;
p=7;

%nondimensional conversion values 
epsilon = eta/delta;
a = phi/(gamma*delta);
b = psi/delta;
c = nu_1/delta;
d = (nu_2*gamma)/(mu*eta);
f = xi * c;

%Initial conditions
IC_N = .2;
IC_A = .03;
IC_E = .8;
IC_hill = [IC_N IC_A IC_E];
%calculating My full model

%Solving NAE-model using ode23
[IVsol_hill, DVsol_hill] = ode45(@(t, y) DEdef_hill(t, y, a,b,c,f,d,p,epsilon), domain, IC_hill);
N_sol_hill = DVsol_hill(:, 1)*gamma;
A_sol_hill = DVsol_hill(:, 2)*gamma;
E_sol_hill = DVsol_hill(:, 3)*mu;



%Solution plot of Hill-model
 fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_hill, N_sol_hill, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_hill) * 1.2]);
ylabel('nutrients','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_hill, A_sol_hill, 'color', algaecolordet, 'linewidth', 3);
hold on;
plot(IVsol_hill, E_sol_hill, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([A_sol_hill; E_sol_hill]) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('algae & EPS','FontSize',20,'Color','k');

% Set common properties
xlim([0, n]);
xlabel('time (days)','FontSize',20,'Color','k');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors


% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)




%Defining NAE-model with hill function inflow/outflow of nutrient rates
function [Dode] = DEdef_hill(I,D,a,b,c,f,d,p,epsilon)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = ((a)/(1+E)^p - (c*N*A)/(1 + N) - (b*N)/(1+E)^p)/epsilon;
dAdt = ((f*N*A)/(1 + N) - A)/epsilon;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
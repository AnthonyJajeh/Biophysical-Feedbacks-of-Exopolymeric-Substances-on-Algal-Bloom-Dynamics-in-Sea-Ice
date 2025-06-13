%Author:Anthony Jajeh
%Date: Febuary 5th, 2025
% Steady state analysis of NAE-model using exponential function as
% inflow/outflow of nutrients functions
clear all; clc; close all;
n=50;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS
% 
% % %Parameter values for fig 4a
% phi = .0001;
% psi = .1;
% mu = .08;
% gamma = .01; 
% nu_1 = .2; 
% nu_2 = .5; 
% xi = .2;
% delta = .007; 
% eta = .03;


% %Parameter values for fig 4b
% phi = .1;
% psi = .5;
% mu = .08;
% gamma = .01; 
% nu_1 = .2; 
% nu_2 = .05; 
% xi = .2;
% delta = .007; 
% eta = .03;

%Parameter values for fig 4c
phi = .01;
psi = .001;
mu = .08;
gamma = .01; 
nu_1 = .2; 
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;

%nondimensional conversion values 
epsilon = eta/delta;
a = phi/(gamma*delta);
b = psi/delta;
c = nu_1/delta;
d = (nu_2*gamma)/(mu*eta);
f = xi * c;

%Initial conditions
IC_N = 5;
IC_A = .03;
IC_E = .01;

IC_exp = [IC_N IC_A IC_E];
value = (IC_A/epsilon)*((f*IC_N)/(IC_N+1)-1);
fprintf('N = %.2f mg N/L\n', value);
%calculating NAE-model solution plots 
[IVsol_exp, DVsol_exp] = ode45(@(t, y) DEdef_exp(t, y, a,b,c,f,d,epsilon), domain, IC_exp);
N_sol_exp = DVsol_exp(:, 1)*gamma;
A_sol_exp = DVsol_exp(:, 2)*gamma;
E_sol_exp = DVsol_exp(:, 3)*mu;



%plotting solution curves of full-model 
% Create a new figure
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_exp, N_sol_exp, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_exp) * 1.2]);
ylabel('nutrients (mg N/L)','FontSize',20,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_exp, A_sol_exp, 'color', algaecolordet, 'linewidth', 3);
hold on;
plot(IVsol_exp, E_sol_exp, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([A_sol_exp; E_sol_exp]) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('algae (mg chl A/L) & EPS (mg XGEQUIV/L)','FontSize',20,'Color','k');

% Set common properties
xlim([0, n]);
xlabel('time (days)','FontSize',20,'Color','k');
set(gca, 'fontsize', 17, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors


% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

%Defining NAE-model
function [Dode] = DEdef_exp(I,D,a,b,c,f,d,epsilon)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = ((a*exp(-E)-(c*A*N)/(N+1)-b*N*exp(-E)))/epsilon;
dAdt = ((f*N*A)/(1 + N) - A)/epsilon;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
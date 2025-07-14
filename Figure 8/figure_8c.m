%Author:Anthony Jajeh
%Date: Febuary 12th, 2025
%Testing to see when algal and EPS blooms appear for various values of
%parameter c
clear all; clc; close all;
n=250;
n_amount = 20;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS
%Parameter values 
phi = .001;
psi = .001;
mu =  .0008;
gamma = .01; 
nu_1_vec = linspace(.01,.4,20);
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;


%nondimensional conversion values 
epsilon = eta/delta;
a = phi/(gamma*delta);
b = psi/delta;
c_vec = nu_1_vec/delta;
d = (nu_2*gamma)/(mu*eta);
f = 5.7;

%Initial conditions
IC_N = 5;
IC_A = .03;
IC_E = .01;


%Allocting space for the maximum values of algae, nutrients, and EPS 
A_max = zeros(1, n_amount);
N_max = zeros(1,n_amount);
E_max = zeros(1,n_amount);

%runs a solution plot for different values of specified parameter 
for i = 1:n_amount
    c = c_vec(i);
    IC_exp = [IC_N IC_A IC_E];
    % Solve model for each c 
    [IVsol_exp, DVsol_exp] = ode23(@(t, y) DEdef_exp(t, y, a,b,c,f,d,epsilon), domain, IC_exp);
    N_sol_exp = DVsol_exp(:, 1)*gamma;
    A_sol_exp = DVsol_exp(:, 2)*gamma;
    E_sol_exp = DVsol_exp(:, 3)*mu;
    
    %Max values of each state variable
    N_max(i)=max(N_sol_exp);
    A_max(i) = max(A_sol_exp);
    E_max(i) = max(E_sol_exp);
    
  
    %plotting solution curves of NAE-model 
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
title("c=",c)

% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

end


fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

yyaxis left;
plot(c_vec, N_max, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_max) * 1.2]);
ylabel('max nutrients (mg N/L)','FontSize',17,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(c_vec, A_max, 'color', algaecolordet, 'linewidth', 3);
hold on;
plot(c_vec, E_max, 'color', EPScolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([max(A_max); max(E_max)]) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('max algae (mg chl-a/L) & EPS (mg XG/L)','FontSize',17,'Color','k');

xlabel('\it{c} (rate that algae uptake nutrients)', 'FontSize', 17);
xlim([min(c_vec),max(c_vec)])
% Add legend
set(gca, 'YColor', 'k'); % <-- Apply black color to right y-axis
legend('Max Nutrients', 'Max Algae', 'Max EPS', 'Location', 'northeast');
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
dNdt =(a*exp(-E)-(c*A*N)/(N+1)-b*N*exp(-E))/epsilon;
dAdt = ((f*N*A)/(1 + N) - A)/epsilon;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
%Author:Anthony Jajeh
%Date: Febuary 12th, 2025
%Testing to see when algal and EPS blooms appear for various values of N_0 
clear all; clc; close all;
n=250;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS

a = 8; %infow of nutrients
b = .5; %outflow of nutrients
c = .8; %Nutrient uptake by algae 
c_p = 1.1; %algal growth rate
d = .5; %EPS growth rate due to algae 


IC_N_values = linspace(0, 100, 20); %Range of initial conditions of nutrients
IC_A = 1; %Initial condition of algae
IC_E = 1; %Initial condition of EPS

%Allocting space for the maximum values of algae, nutrients, and EPS 
A_max = zeros(1,length(IC_N_values));
N_max = zeros(1,length(IC_N_values));
E_max = zeros(1,length(IC_N_values));

%runs a solution plot for different values of specified parameter 
for i = 1:length(IC_N_values)
    IC_N = IC_N_values(i);
    
    IC_exp = [IC_N IC_A IC_E];
    % Solve simplified model for current N
    
    [IVsol_exp, DVsol_exp] = ode23(@(t, y) DEdef_exp(t, y, a,b,c,c_p,d), domain, IC_exp);
    N_sol_exp = DVsol_exp(:, 1);
    A_sol_exp = DVsol_exp(:, 2);
    E_sol_exp = DVsol_exp(:, 3);
    
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


% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

end

figure;
% Plot N_max vs b_values
plot(IC_N_values, N_max, 'o-', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'N_{max}')
xlabel('N_0 values')
ylabel('Maximum values of Nutrients')
title('Maximum N vs N_0 values')
grid on
figure
% Plot A_max vs b_values
plot(IC_N_values, A_max, 's-', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'A_{max}')
xlabel('N_0 values')
ylabel('Maximum values of Algae')
title('Maximum A vs N_0 values')
grid on
figure
% Plot A_max vs b_values
plot(IC_N_values, E_max, 's-', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'E_{max}')
xlabel('N_0 values')
ylabel('Maximum values of EPS')
title('Maximum E vs N_0 values')
grid on

%Defining NAE-model
function [Dode] = DEdef_exp(I,D,a,b,c,c_p,d)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt =a*exp(-E)-(c*A*N)/(N+1)-b*N*exp(-E);
dAdt = (c_p*N*A)/(1 + N) - A;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
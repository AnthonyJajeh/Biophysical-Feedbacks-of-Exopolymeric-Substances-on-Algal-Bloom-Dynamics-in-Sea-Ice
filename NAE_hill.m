%Author:Anthony Jajeh
%Date: Febuary 5th, 2025
% Steady state analysis of NAE-model using hill function as
% inflow/outflow of nutrients functions
clear all; clc; close all;
n=250;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 
%Parameter values 
a = 8; %infow of nutrients
b = .1; %outflow of nutrients
c = .8; %Nutrient uptake by algae 
c_p = 1.3; %algal growth rate
d = .5; %EPS growth rate due to algae 
p=1; %negative feedback

%creating vectors
IC_N = 15;
IC_A = 1;
IC_E = 1;
IC_hill = [IC_N IC_A IC_E];
%calculating My full model

%Solving NAE-model using ode23
[IVsol_hill, DVsol_hill] = ode23(@(t, y) DEdef_hill(t, y, a,b,c,c_p,d,p), domain, IC_hill);
N_sol_hill = DVsol_hill(:, 1);
A_sol_hill = DVsol_hill(:, 2);
E_sol_hill = DVsol_hill(:, 3);

%Creating vectors for storage for a and p values for a-p bifurcation plot
Avec = linspace(.01,10,n);
avec = zeros(0,n);
pvec = zeros(0,n);

%Implicitely solving for p using realchar function defined below then
%updating a parameter as a function of algae
for i=1:n
    A_sim = Avec(i);
    pvec(i) = realchar(A_sim,b,c,c_p,d,10);
    % for each value of A we implicitly solve for p
    avec(i) = (A_sim*c*(c_p - 1)*(A_sim*d + 1)^p + b*c_p)/((c_p - 1)*c_p);
end
%two parameter a-p bifurcation plot
 figure;
 plot(avec,pvec,'LineWidth',2,'color','k')
 xlabel('a (inflow of nutrients)','FontSize',25);
 ylabel('p (negative feedback)','FontSize',25);
 text(4, 1, 'Spiral sink', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
 text(4, 6, 'Periodic ', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
 ylim([0,15])
 xlim([0,15])
%Solution plot of NAE-model
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




%implicitly solve for p
function x = realchar(A,b,c,c_p,d,x0)
re = @(p)(3*(c_p)^2*(A*d + 1)^2*b*((c*A + 1/3)*(c_p)^2 - 2*c*A*(c_p) + c*A)*(A*d + 1)^(-1 - p) + b^2*(c_p)^4*(A*d + 1)^2*(A*d + 1)^(-1 - 2*p) + 2*A*((1/2 + (A)^2*d*c + (c + (-p/2 + 1/2)*d)*A)*(c_p)^2 - 2*A*c*(A*d + 1)*c_p + A*c*(A*d + 1))*c*(c_p - 1)^2)/((c_p)^4*(A*d + 1));
x = fzero(re,x0);
end

%Defining NAE-model with hill function inflow/outflow of nutrient rates
function [Dode] = DEdef_hill(I,D,a,b,c,c_p,d,p)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = (a)/(1+E)^p - (c*N*A)/(1 + N) - (b)/(1+E)^p;
dAdt = (c_p*N*A)/(1 + N) - A;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
%Author:Anthony Jajeh 
%Date: Febuary 5th, 2024
%slow-Quasi-steady-state figures by using MATLAB built
%Implicitly solving for two parameter bifurcation regime of c_p and b to
%determine the stability of the equilibria
%in ODE solver 

clear all; clc; close all;

%Parameters
n=10;
domain = [0 n];
a = 3;
c = 1.1;
d = 1.3;

EPScolordet = 1/255*[125,91,166]; % color for EPS graph


% Solving for b as a function of a,c,c_p, and d to determine the stability of equilibrium point 
% If the inside of the LambertW= -1/e then LambertW(-1/e)=-1
%Inside function comes from inside the LambertW function of the equilibrium
%point

syms c_p b
inside = (c_p * d * (a * c_p - a - b)) / ((c * (c_p - 1)))+1/exp(1);

b_solution = solve(inside, b);
% Convert symbolic solution to a function for plotting
b_func = matlabFunction(b_solution, 'Vars', c_p);

% Plot b as a function of c_p
c_p_vals = linspace(1, 10, 500);
b_vals = b_func(c_p_vals);
%Two parameter bifurcation graph of c_p vs b
figure;
plot(c_p_vals, b_vals, 'LineWidth', 2,'Color','k');
xlabel('$\overline{c}$ (growth rate of algae)','interpreter','latex','FontSize', 25);
ylabel('b (outflow of nutrients)','FontSize',25);
text(3, 15, 'Unstable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
text(6, 5, 'Stable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
ylim([0,max(b_vals)])
xlim([1,max(c_p_vals)])

%Solving dEdt using ode15s
b= 2;
c_p = 4;
IC_E = 4 ; 
[IVsol_slow, DVsol_slow] = ode15s(@(t, E) DEdef(t, E, a,b,c,c_p,d), domain, IC_E);

% Solution plot for dEdt
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

plot(IVsol_slow, DVsol_slow, 'color', EPScolordet, 'linewidth', 3);
ylim([0, max(DVsol_slow) * 1.2]);
ylabel('EPS','Color','k');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Set common properties
xlim([0, n]);
xlabel('time (days)');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Add legend
legend('EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)


hold off;
% Set the figure size and save as PDF, PNG, and FIG
set(gcf, 'units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]); % Set the paper to have width 6 and height 4

%Defining the dEdt ODE
function [dEdt] = DEdef(~,E,a,b,c,c_p,d)

%set of odes
dEdt = (d*c_p*(c_p*a-a-b))/(c*(c_p-1))*exp(-E)-E;
end

%Author:Anthony Jajeh 
%Date: June 5th, 2024
%slow-Quasi-steady-state figures by using MATLAB built
%Implicitly solving for two parameter bifurcation regime of c_p and b to
%determine the stability of the equilibria
%in ODE solver 

clear all; clc; close all;
EPScolordet = 1/255*[125,91,166]; % Color for EPS graph

%Parameters
n=25;
domain = [0 n];

%Parameter values for fig 2
phi = .001;
psi = .001;
mu = .08;
gamma = .01; 
nu_1 = .2; 
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;


%nondimensional conversion values 
c = nu_1/delta;
d = (nu_2*gamma)/(mu*eta);
f = xi * c;


f = 3; %Inflow of nturients
c = 1.1; %Nutrient uptake by algae
d = 1.3; %Growth rate of EPS due to algae

% Solving for b as a function of a,c,f, and d to determine the stability of equilibrium point 
% If the inside of the LambertW= -1/e then LambertW(-1/e)=-1
%Inside function comes from inside the LambertW function of the equilibrium
%point

syms a b 
inside = (f * d * (a * f - a - b)) / ((c * (f - 1)))+1/exp(1);

b_solution = solve(inside, b);
% Convert symbolic solution to a function for plotting
b_func = matlabFunction(b_solution, 'Vars', a);

% Plot b as a function of a
a_vals = linspace(1, 10, 500);
b_vals = b_func(a_vals);
%Two parameter bifurcation graph of a vs b
figure;
plot(a_vals, b_vals, 'LineWidth', 2,'Color','k');
xlabel('\it{a} (growth rate of algae)','FontSize', 25);
ylabel('\it{b} (outflow of nutrients)','FontSize',25);
text(3, 15, 'Unstable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
text(6, 5, 'Stable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
ylim([0,max(b_vals)])
xlim([1,max(a_vals)])

%Solving dEdt using ode15s
%nondimensional conversion values 
a = phi/(gamma*delta);
b = psi/delta;
IC_E = .79 ; 
[IVsol_slow, DVsol_slow] = ode15s(@(t, E) DEdef(t, E, a,b,c,f,d), domain, IC_E);

% Solution plot for dEdt
fig = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

plot(IVsol_slow, DVsol_slow* mu, 'color', EPScolordet, 'linewidth', 3);
ylim([0, max(DVsol_slow*mu) * 1.2]);
ylabel('EPS (mg XG/L)','Color','k');
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
function [dEdt] = DEdef(~,E,a,b,c,f,d)

%set of odes
dEdt = (d*f*(f*a-a-b))/(c*(f-1))*exp(-E)-E;
end

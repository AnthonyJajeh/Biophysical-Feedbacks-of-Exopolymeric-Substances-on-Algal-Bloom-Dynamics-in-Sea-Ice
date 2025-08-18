%Author:Anthony Jajeh 
%Date: Febuary 5th, 2024
%Quasi-steady-state figures by using MATLAB built
%in ODE solver
clear all;clc;close all;

%domain 
n=500;
domain = [0 n];

% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 


%Parameter values for fig 3
phi = .01;
psi = .01;
mu = .000008;
gamma = .01; 
nu_1 = .2; 
nu_2 = .05; 
xi = .2;
delta = .007; 
eta = .03;
sigma = .8; % example scaling factor, adjust as needed


%nondimensional conversion values 
epsilon = eta/delta;
a = phi/(gamma*delta);
b = psi/delta;
c = nu_1/delta;
d = (nu_2*gamma)/(mu*eta);
f = xi * c;
h = (sigma*gamma)/mu;

%Initial conditions
IC_N = 10;
IC_A = .03;
IC_E = .01;

%initial condition vector 
IC_fast = [IC_N IC_A];
IC_slow = IC_E;
IC_full = [IC_N IC_A IC_E];

%Solving QSS fast-model 
[IVsol_fast, DVsol_fast] = ode23s(@(t, y) DEdef_fast(t, y, a,b,c,f,IC_E), domain, IC_fast);
N_sol_fast = DVsol_fast(:, 1)*gamma;
A_sol_fast = DVsol_fast(:, 2)*gamma;

%Solving QSS slow-model
[IVsol_slow, DVsol_slow] = ode45(@(t, E) DEdef_slow(t, E, a,b,c,f,d), domain, IC_E);


%calculating full-model solution plots 
[IVsol_exp, DVsol_exp] = ode23s(@(t, y) DEdef_exp(t, y, a,b,c,f,d,epsilon), domain, IC_full);
N_sol_exp = DVsol_exp(:, 1)*gamma;
A_sol_exp = DVsol_exp(:, 2)*gamma;
E_sol_exp = DVsol_exp(:, 3)*mu;

%Solving tracking-model using ode23
[IVsol_tracking, DVsol_tracking] = ode23s(@(t, y) DEdef_tracking(t, y, a,b,c,f,h), domain, IC_fast);
N_sol_tracking = DVsol_tracking(:, 1)*gamma;
A_sol_tracking = DVsol_tracking(:, 2)*gamma;

E_plot = IC_E * ones(size(E_sol_exp));

fig = figure;
hold on;

% Plot EPS from slow model
plot(IVsol_slow, DVsol_slow*mu, 'Color', EPScolordet, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Slow Model');

% Plot EPS from full model
plot(IVsol_exp, E_sol_exp, 'Color', EPScolordet, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Full Model');

%Plot EPS from fast-model
plot(IVsol_exp, E_plot*mu, 'Color', EPScolordet, 'LineWidth', 2,'LineStyle', ":", 'DisplayName', 'Fast Model');

sigma = .8; % example scaling factor, adjust as needed
E_tracking = sigma * DVsol_tracking(:, 2)*gamma;
%Plot EPS from tracking model E(t) = sigma * A(t)
plot(IVsol_tracking, E_tracking, 'Color', EPScolordet, 'LineWidth', 2,'LineStyle', '-.', 'DisplayName', 'Tracking Model');

% Set axes limits dynamically based on max EPS value from all
maxnutrient = max([DVsol_slow*mu; E_sol_exp; E_tracking]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('Time (days)', 'Color', 'k');
ylabel('EPS (mg XG/L)', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);

fig = figure;
hold on;

% Plot algae from full model
plot(IVsol_exp, A_sol_exp, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', '-',  'DisplayName', 'Full Model');

% Plot algae from fast model
plot(IVsol_fast, A_sol_fast, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', ':',  'DisplayName', 'Fast Model');

% Plot aglae from Tracking model 
plot(IVsol_tracking, A_sol_tracking, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', '-.',  'DisplayName', 'Tracking Model');

% Set axes limits dynamically based on max algae value from all
maxnutrient = max([A_sol_exp; A_sol_tracking; A_sol_fast]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('Time (days)', 'Color', 'k');
ylabel('Algae (mg chl-a/L)', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);


%Plot for figure 3c
fig = figure;
hold on;

% Plot nutrient from full model
plot(IVsol_exp, N_sol_exp, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', '-',  'DisplayName', 'Full Model');

% Plot nutrient from fast model
plot(IVsol_fast, N_sol_fast, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', ':',  'DisplayName', 'Fast Model');

% Plot nutrient from Tracking model 
plot(IVsol_tracking, N_sol_tracking, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', '-.',  'DisplayName', 'Tracking Model');

% Set axes limits dynamically based on max algae value from all
maxnutrient = max([N_sol_exp; N_sol_tracking; N_sol_fast]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('Time (days)', 'Color', 'k');
ylabel('Nutrients (mg N/L)', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);







%Defining QSS-fast-model
function [Dode] = DEdef_fast(I,D,a,b,c,f,E_0)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);

%set of odes
dNdt = a*exp(-E_0)-(c*N*A)/(N+1) - b*N*exp(-E_0);
dAdt = (f*N*A)/(N+1)-A;

% odes in vector form
Dode = [dNdt; dAdt];
end

%Defining the dEdt ODE
function [dEdt] = DEdef_slow(~,E,a,b,c,f,d)

%set of odes
dEdt = (d*f*(f*a-a-b))/(c*(f-1))*exp(-E)-E;
end


%Defining NAE-full
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

%Defining full model with E(t) = sigma * A(t)
function [Dode] = DEdef_tracking(I,D,a,b,c,f,h)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);


%set of odes
dNdt = a*exp(-h*A)-(c*N*A)/(1+N)-b*N*exp(-h*A);
dAdt = (f*N*A)/(1 + N) - A;

% odes in vector form
Dode = [dNdt; dAdt];
end
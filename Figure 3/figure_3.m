%Author:Anthony Jajeh 
%Date: Aug 18, 2025
%Quasi-steady-state figures by using MATLAB built
%in ODE solver
clear all;clc;close all;

%domain 
n=300;
domain = [0 n];

% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 


% %Parameter values for fig 3 a/b<
% phi = .0001;
% psi = .5;
% mu = .001;
% gamma = .01; 
% nu = .2; 
% rho = .75; 
% xi = .2;
% delta = .007; 
% eta = .03;
% sigma = .8; % example scaling factor, adjust as needed

% % %Parameter values for fig 3 a/b>
phi = .0001;
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;
sigma = .8; % example scaling factor, adjust as needed

%nondimensional conversion values 
epsilon = eta/delta;
a = phi/(gamma*delta);
b = psi/delta;
c = nu/delta;
d = (rho*gamma)/(mu*eta);
f = xi * c;
h = (sigma*gamma)/mu;

%Initial conditions
IC_N = .005/gamma;
IC_A = .003/gamma;
IC_E = .001/mu;

%initial condition vector 
IC_fast = [IC_N IC_A];
IC_slow = IC_E;
IC_full = [IC_N IC_A IC_E];

opts_slow = odeset('RelTol',1e-7,'AbsTol',1e-9,'InitialStep',1e-4,'MaxStep',1, ...
                   'NonNegative',[1 1 1]);
%Solving QSS fast-model 
[IVsol_slow, DVsol_slow] = ode15s(@(t, y) DEdef_slow(t, y, a,b,c,f,d,IC_E), domain, IC_full,opts_slow);
N_sol_slow = DVsol_slow(:, 1);
A_sol_slow = DVsol_slow(:, 2);
E_sol_slow = DVsol_slow(:,3);


%calculating full-model solution plots 
[IVsol_exp, DVsol_exp] = ode23s(@(t, y) DEdef_exp(t, y, a,b,c,f,d,epsilon), domain, IC_full);
N_sol_exp = DVsol_exp(:, 1);
A_sol_exp = DVsol_exp(:, 2);
E_sol_exp = DVsol_exp(:, 3);

%Solving fast-model using ode23
[IVsol_fast, DVsol_fast] = ode23s(@(t, y) DEdef_fast(t, y, a,b,c,f,h), domain, IC_fast);
N_sol_fast = DVsol_fast(:, 1);
A_sol_fast = DVsol_fast(:, 2);

E_plot = IC_E * ones(size(E_sol_exp));

fig1 = figure;
hold on;

% Plot EPS from slow model
plot(IVsol_slow, DVsol_slow(:,3), 'Color', EPScolordet, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Regime 1');

sigma = .8; % example scaling factor, adjust as needed
E_tracking = sigma * DVsol_fast(:, 2);
%Plot EPS from tracking model E(t) = sigma * A(t)
plot(IVsol_fast, E_tracking, 'Color', EPScolordet, 'LineWidth', 2,'LineStyle', '-.', 'DisplayName', 'Regime 2');


% Plot EPS from full model
plot(IVsol_exp, E_sol_exp, 'Color', EPScolordet, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Regime 3');


% Set axes limits dynamically based on max EPS value from all
maxnutrient = max([E_sol_exp; DVsol_slow(:,3)]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('time', 'Color', 'k');
ylabel('EPS', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);

fig2 = figure;
hold on;


% Plot aglae from slow model 
plot(IVsol_fast, A_sol_fast, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', '-.',  'DisplayName', 'Regime 1');

% Plot algae from fast model
plot(IVsol_slow, A_sol_slow, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', ':',  'DisplayName', 'Regime 2');

% Plot algae from full model
plot(IVsol_exp, A_sol_exp, 'Color', algaecolordet, 'LineWidth', 2,'LineStyle', '-',  'DisplayName', 'Regime 3');


% Set axes limits dynamically based on max algae value from all
maxnutrient = max([A_sol_exp; A_sol_fast; A_sol_slow]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('time', 'Color', 'k');
ylabel('algae', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);


%Plot for figure 3c
fig3 = figure;
hold on;


% Plot nutrient from fast model
plot(IVsol_slow, N_sol_slow, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', ':',  'DisplayName', 'Regime 1');

% Plot nutrient from Tracking model 
plot(IVsol_fast, N_sol_fast, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', '-.',  'DisplayName', 'Regime 2');

% Plot nutrient from full model
plot(IVsol_exp, N_sol_exp, 'Color', nutrientcolordet, 'LineWidth', 2,'LineStyle', '-',  'DisplayName', 'Regime 3');

% Set axes limits dynamically based on max algae value from all
maxnutrient = max([N_sol_exp; N_sol_fast; N_sol_slow]);
ylim([0, maxnutrient * 1.2]);
xlim([0, n]);

xlabel('time', 'Color', 'k');
ylabel('nutrients', 'Color', 'k');
set(gca, 'FontSize', 20, 'XColor', 'k', 'YColor', 'k');

legend('Location', 'northeast');
legend boxoff;

hold off;

set(gcf, 'Units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]);






fname1 = 'fig3f';
fname2= 'fig3e';
fname3 = 'fig3d';
nice_graphing(fname1,fig1)
nice_graphing(fname2, fig2)
nice_graphing(fname3,fig3)

function nice_graphing(fname, fig)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .8; % feel free to play with this ratio
set(findall(fig,'-property','FontSize'),'FontSize',24) % adjust fontsize to your document
set(findall(fig,'-property','Box'),'Box','on') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
lgd = findall(fig, 'Type', 'Legend');
    set(lgd, 'Box', 'off');     % ensure no border if a legend exists
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(fig, strcat(fname,'.png'), 'ContentType', 'vector');
end


%Defining QSS-slow-model
function [Dode] = DEdef_slow(I,D,a,b,c,f,d,E_0)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = a*exp(-E_0)-(c*N*A)/(N+1) - b*N*exp(-E_0);
dAdt = (f*N*A)/(N+1)-A;
dEdt = (d*f*(a*f-a-b))/(exp(E)*c*(f-1));
% odes in vector form
Dode = [dNdt; dAdt;dEdt];
end

%Defining full-model
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
function [Dode] = DEdef_fast(I,D,a,b,c,f,h)
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
%Author:Anthony Jajeh 
%Date: Aug 18, 2025
%Figure 1 by using MATLAB built
%in ODE solver
clear all;clc;close all;



% Define colors for deterministic results
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 


% %Parameter values for fig 1ac trivial
phi_0 = .0001;
psi_0 = .05;
mu_0= .001;
gamma_0 = .01; 
nu_0 = .2; 
rho_0 = .75; 
xi_0 = .2;
delta_0 = .007; 
eta_0 = .03;% 

%Parameter values for fig 1bc nontrivial
 phi_1 = .01;
 psi_1 = .01;
 mu_1 = .001;
 gamma_1 = .01; 
 nu_1 = .2; 
 rho_1 = .75; 
 xi_1 = .2;
 delta_1 = .007; 
 eta_1 = .03;


%nondimensional conversion values 
a_0 = phi_0/(gamma_0*delta_0);
b_0 = psi_0/delta_0;
c_0 = nu_0/delta_0;
d_0 = (rho_0*gamma_0)/(mu_0*eta_0);
f_0 = xi_0 * c_0;
epsilon_0 = eta_0/delta_0;


%nondimensional conversion values 
a_1 = phi_1/(gamma_1*delta_1);
b_1 = psi_1/delta_1;
c_1 = nu_1/delta_1;
d_1 = (rho_1*gamma_1)/(mu_1*eta_1);
f_1 = xi_1 * c_1;
epsilon_1 = eta_1/delta_1;


%domain 
n=10;
domain = [0 n];

%Initial conditions
IC_N = .2/gamma_0;
IC_A = .0002/gamma_0;
IC_E = .002/mu_0;



%initial condition vector 
IC_fast = [IC_N IC_A];
IC_slow = IC_E;
%scaled values

%Solving QSS fast-model 
[IVsol_fast_triv, DVsol_fast_triv] = ode23s(@(t, y) DEdef_fast(t, y, a_0,b_0,c_0,f_0,IC_E), domain, IC_fast);
N_sol_fast_triv = DVsol_fast_triv(:, 1);
A_sol_fast_triv = DVsol_fast_triv(:, 2);

opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-2, 'NonNegative',1);
[IVsol_slow_triv, DVsol_slow_triv] = ode15s(@(t,y) DEdef_E(t,y), domain, IC_E, opts);

%Solving QSS fast-model 
[IVsol_fast_non, DVsol_fast_non] = ode23s(@(t, y) DEdef_fast(t, y, a_1,b_1,c_1,f_1,IC_E), domain, IC_fast);
N_sol_fast_non = DVsol_fast_non(:, 1);
A_sol_fast_non = DVsol_fast_non(:, 2);

opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-2, 'NonNegative',1);
[IVsol_slow_non, DVsol_slow_non] = ode15s(@(t,y) DEdef_slow(t,y,a_1,b_1,c_1,f_1,d_1), domain, IC_E, opts);


NP_0=a_0/b_0;
AP_0=0;
NP_1=1/(f_1-1);
AP_1= f_1*(a_1*f_1-a_1-b_1)/(exp(IC_E)*c_1*(f_1-1));
%Create a new figure

%Trivial
% Create a new figure
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figsfast_triv = figure;
set(figsfast_triv, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis 
yyaxis left;
plot(IVsol_fast_triv, N_sol_fast_triv, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_fast_triv)*1.2]);
ylabel('nutrients','FontSize',20,'Color','k','Interpreter','latex');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_fast_triv, A_sol_fast_triv, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_fast_triv, IC_E*ones(size(IVsol_fast_triv)), 'color', EPScolordet,'LineStyle', '--','linewidth', 3);
ylim([0, max(max(IC_E), max(A_sol_fast_triv))*1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');

% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
%Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figslow_triv = figure;
set(figslow_triv, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_slow_triv, NP_0*ones(size(IVsol_slow_triv)), 'color', nutrientcolordet, 'LineStyle', '--', 'linewidth', 3);
ylim([0, max(NP_0)*1.2]);
ylabel('nutrients','FontSize',20,'Color','k','Interpreter','latex');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
% Your original plot order (keep it)
h_alg = plot(IVsol_slow_triv, AP_0*ones(size(IVsol_slow_triv)), ...
    'Color', algaecolordet, 'LineStyle', '--', 'LineWidth', 3);
hold on
h_eps = plot(IVsol_slow_triv, DVsol_slow_triv, ...
    'Color', EPScolordet, 'LineStyle', '-', 'LineWidth', 3);

% Re-plot algae ON TOP, but hide from legend
h_alg_top = plot(IVsol_slow_triv, AP_0*ones(size(IVsol_slow_triv)), ...
    'Color', algaecolordet, 'LineStyle', '--', 'LineWidth', 3);

% Hide the top copy from legend (two robust ways; either is fine)
h_alg_top.HandleVisibility = 'off';
% OR:
% h_alg_top.Annotation.LegendInformation.IconDisplayStyle = 'off';

ylim([0, max(DVsol_slow_triv)*1.2]);
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');

hold off
% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
%Add legend
legend('nutrients', 'algae','EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)=
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors


%Trivial
% Create a new figure
figsfast_non = figure;
set(figsfast_non, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis 
yyaxis left;
plot(IVsol_fast_non, N_sol_fast_non, 'color', nutrientcolordet, 'linewidth', 3);
ylim([0, max(N_sol_fast_non)*1.2]);
ylabel('nutrients','FontSize',20,'Color','k','Interpreter','latex');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_fast_non, A_sol_fast_non, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
plot(IVsol_fast_non, IC_E*ones(size(IVsol_fast_non)), 'color', EPScolordet,'LineStyle', '--','linewidth', 3);
ylim([0, max(max(IC_E), max(A_sol_fast_non))*1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');


% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
%Add legend
%legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors

figslow_non = figure;
set(figslow_non, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(IVsol_slow_non, NP_1*ones(size(IVsol_slow_non)), 'color', nutrientcolordet, 'LineStyle', '--', 'linewidth', 3);
ylim([0, max(IVsol_slow_non)*1.2]);
ylabel('nutrients','FontSize',20,'Color','k','Interpreter','latex');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
plot(IVsol_slow_non, AP_1*ones(size(IVsol_slow_non)), 'color', algaecolordet, 'LineStyle', '--', 'linewidth', 3);
plot(IVsol_slow_non, DVsol_slow_non, 'color', EPScolordet,'LineStyle', '-', 'linewidth', 3);
ylim([0, max(max(IC_E), max(A_sol_fast_non))*1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('$\mathrm{algae\ \&\ EPS}$','FontSize',20,'Color','k','Interpreter','latex');


% Set common properties
xlim([0, n]);
xlabel('time','FontSize',20,'Color','k');
% %Add legend
%legend('nutrients', 'algae','EPS', 'Location', 'northeast');
hold off
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis text and tick colors

fname1 = 'fig1a';
fname2 = 'fig1b';
fname3= 'fig1c';
fname4 = 'fig1d';
nice_graphing(fname1,figsfast_triv)
nice_graphing(fname2,figsfast_non)
nice_graphing(fname3,figslow_triv)
nice_graphing(fname4,figslow_non)




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
exportgraphics(fig, strcat(fname,'.png'), 'Resolution', 300);
exportgraphics(fig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(fig, strcat(fname,'.fig'));
end
% 
% exportgraphics(figsfast_triv, strcat(fname1,'.png'), 'Resolution', 300);
% exportgraphics(figsfast_triv, strcat(fname1,'.pdf'), 'ContentType', 'vector');
% saveas(figsfast_triv, strcat(fname1,'.fig'));
% 
% exportgraphics(figsfast_non, strcat(fname2,'.png'), 'Resolution', 300);
% exportgraphics(figsfast_non, strcat(fname2,'.pdf'), 'ContentType', 'vector');
% saveas(figsfast_non, strcat(fname2,'.fig'));
% 
% exportgraphics(figslow_triv, strcat(fname3,'.png'), 'Resolution', 300);
% exportgraphics(figslow_triv, strcat(fname3,'.pdf'), 'ContentType', 'vector');
% saveas(figslow_triv, strcat(fname3,'.fig'));
% 
% exportgraphics(figslow_non, strcat(fname4,'.png'), 'Resolution', 300);
% exportgraphics(figslow_non, strcat(fname4,'.pdf'), 'ContentType', 'vector');
% saveas(figslow_non, strcat(fname4,'.fig'));

%Defining QSS-NA-model
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
dEdt = (d*f*(a*f-a-b))/(exp(E)*c*(f-1))-E;
end
% %THIS IS SUPPOSEED TO BE EXPONENTIAL DECAY FOR TRIVIAL SOLN
% function [dEdt] = DEdef_slow(~,E,a,b,c,f,d)
% 
% %set of odes
% dEdt = (d*f*(a*f-a-b))/(exp(E)*c*(f-1))-E;
% end

%Defining the dEdt ODE
function [dEdt] = DEdef_E(~,E)

%set of odes
dEdt = -E;
end
% %THIS IS SUPPOSEED TO BE EXPONENTIAL DECAY FOR TRIVIAL SOLN
% function [dEdt] = DEdef_slow(~,E,a,b,c,f,d)
% 
% %set of odes
% dEdt = (d*f*(a*f-a-b))/(exp(E)*c*(f-1))-E;
% end

%Author:Anthony Jajeh
%Date: Aug 26, 2025
%Testing to see when algal and EPS blooms appear for various values of parameter a 
clear all; clc; close all;
n=3000;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS

%Parameter values 
phi_vec = linspace(.001,.1,20);
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;% 


%Initial conditions
IC_N = .005;
IC_A = .003;
IC_E = .001;

%Allocting space for the maximum values of algae, nutrients, and EPS 
A_max = zeros(1,length(phi_vec));
N_max = zeros(1,length(phi_vec));
E_max = zeros(1,length(phi_vec));

%runs a solution plot for different values of specified parameter 
for i = 1:length(phi_vec)
    phi = phi_vec(i);
    IC_exp = [IC_N IC_A IC_E];
    % Solve simplified model for current a
    
    [IVsol_exp, DVsol_exp] = ode23s(@(t, y) DEdef_exp(t, y, phi,psi,nu,xi,delta,rho,eta,mu,gamma), domain, IC_exp);
    N_sol_exp = DVsol_exp(:, 1);
    A_sol_exp = DVsol_exp(:, 2);
    E_sol_exp = DVsol_exp(:, 3);
    
    %Max values of each state variable
    N_max(i) = max(N_sol_exp);
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
title("a=",phi)
% Add legend
legend('Nutrients', 'Algae', 'EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

end


figp = figure;
set(fig, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

% Plot nutrients on the left y-axis
yyaxis left;
plot(phi_vec, N_max, 'color', nutrientcolordet, 'linewidth', 3);
plot(phi_vec, A_max, 'color', algaecolordet, 'linewidth', 3,'LineStyle','-');
ylim([0, max([max(N_max); max(A_max)]) * 1.2]);
ylabel('max nutrients \& algae','FontSize',17,'Color','k');
set(gca, 'YColor', 'k'); % Set the left axis color to black

% Plot algae and EPS on the right y-axis
yyaxis right;
hold on;
plot(phi_vec, E_max, 'color', EPScolordet, 'linewidth', 3);
ylim([0,  max(E_max) * 1.2]); % Ensures that the y-axis accommodates the largest value of algae or EPS
ylabel('max EPS','FontSize',17,'Color','k');

xlabel('$\phi$', 'FontSize', 20);
xlim([min(phi_vec),max(phi_vec)])

set(gca, 'YColor', 'k'); % <-- Apply black color to right y-axis

% Add legend
legend('Max Nutrients', 'Max Algae', 'Max EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)

%Defining NAE-model

fname = 'fig8a';
nice_graphing(fname, figp)

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



%Defining NAE-model
function [Dode] = DEdef_exp(I,D,phi,psi,nu,xi,delta,rho,eta,mu,gamma)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt = phi * exp(-E/mu)- (nu*N*A)/(N+gamma)-psi*N*exp(-E/mu);
dAdt = (xi*nu*N*A)/(N+gamma)-delta*A;
dEdt = rho*A-eta*E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end
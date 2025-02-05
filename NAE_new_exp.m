clear all; clc; close all;
n=250;
domain = [0 n];
algaecolordet = 1/255*[118,176,65]; % color for algae (green)
nutrientcolordet = 1/255*[255,201,20]; % color for nutrients (yellow)\
EPScolordet = 1/255*[125,91,166]; % color for EPS 

a = 10;
b = .1; %Range of b values for plotting
c = .8;
c_p = 1.3;
d = .5;
%creating vectors
IC_N = 15;
IC_A = 1;
IC_E = 1;
IC_exp = [IC_N IC_A IC_E];
%calculating My full model

[IVsol_exp, DVsol_exp] = ode23(@(t, y) DEdef_exp(t, y, a,b,c,c_p,d), domain, IC_exp);
N_sol_exp = DVsol_exp(:, 1);
A_sol_exp = DVsol_exp(:, 2);
E_sol_exp = DVsol_exp(:, 3);

% %creating vectors
 Avec = linspace(.01,5,n);
 avec = zeros(0,n);
 bvec = zeros(0,n);
% 
 for i=1:n
     A_sim = Avec(i);
     avec(i) = -(exp(d*A_sim)*(A_sim*c*(c_p)^2-4*c*A_sim*c_p+3*c*A_sim+(c_p)^2-exp(2*d*A_sim)*sqrt(exp(-4*d*A_sim)*(A_sim^2*(c_p-1)^4*(c)^2+4*(c_p)^2*A_sim*(d*A_sim+.5)*(c_p-1)^2*c+(c_p)^4))))/(2*(c_p-1)*(c_p)^2);
     % for each value of A we implicitly solve for b
     bvec(i) = -(-exp(-d*A_sim)*avec(i)*(c_p)^2+c*A_sim*c_p+avec(i)*exp(-d*A_sim)*c_p-c*A_sim)/(exp(-d*A_sim)*c_p);
 end
% 
 %plotting a_b bif graph
 figure;
 plot(avec,bvec,'LineWidth',2,'color','k')
 xlabel('a (inflow of nutrients)','FontSize',25);
 ylabel('b (outflow of nutrients)','FontSize',25);
 text(10, 1, 'Spiral sink', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
 text(30, .4, 'Periodic ', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
 ylim([0,max(bvec)])
 xlim([0,max(avec)])

%plotting solution curves 
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
%Author:Anthony Jajeh 
%Date: Aug 19th, 2025
%slow-Quasi-steady-state figures by using MATLAB built
%Implicitly solving for two parameter bifurcation regime of c_p and b to
%determine the stability of the equilibria
%in ODE solver 

clear all; clc; close all;
EPScolordet = 1/255*[125,91,166]; % Color for EPS graph

%Parameters
n=10;
domain = [0 n];

%Parameter values for fig 2
phi = .01;
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;



%nondimensional conversion values 
c = nu/delta;
d = (rho*gamma)/(mu*eta);
f = xi * c;

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
fig2 = figure;
plot(a_vals, b_vals, 'LineWidth', 2,'Color','k');
xlabel('\it a \rm (growth rate of algae)','FontSize', 25);
ylabel('\it b \rm (outflow of nutrients)','FontSize',25);
text(1.5, 20, 'Unstable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
text(3, 7, 'Stable', 'Color', 'black', 'FontSize', 20, 'FontWeight', 'bold');
ylim([0,max(b_vals)])
xlim([1,max(a_vals)])

%Solving dEdt using ode15s
%nondimensional conversion values 

phi = .01;
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;

a = phi/(gamma*delta);
b = psi/delta;
c = nu/delta;
d = (rho*gamma)/(mu*eta);
f = xi * c;

IC_E = .001/mu; 
[IVsol_slow, DVsol_slow] = ode23s(@(t, E) DEdef(t, E, a,b,c,f,d), domain, IC_E);

% Solution plot for dEdt
fig1 = figure;
set(fig1, 'defaultAxesColorOrder', [0 0 0; 0 0 0]);
hold on;

plot(IVsol_slow, DVsol_slow, 'color', EPScolordet, 'linewidth', 3);
ylim([0, max(DVsol_slow) * 1.2]);
ylabel('EPS','Color','k');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Set common properties
xlim([0, n]);
xlabel('time');
set(gca, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k'); % Set axis tick label colors to black


% Add legend
legend('EPS', 'Location', 'northeast');
legend boxoff; % Hide the legend's axes (border and background)


hold off;
% Set the figure size and save as PDF, PNG, and FIG
set(gcf, 'units', 'inches', 'Position', [2, 2, 6, 4]);
set(gcf, 'PaperSize', [10, 6]); % Set the paper to have width 6 and height 4

fname1 = 'fig2a';
fname2= 'fig2b';
nice_graphing(fname2, fig2)
nice_graphing(fname1,fig1)

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
exportgraphics(fig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(fig,strcat(fname,'.fig'))

end


%Defining the dEdt ODE
function [dEdt] = DEdef(~,E,a,b,c,f,d)

%set of odes
dEdt = (d*f*(a*f-a-b))/(exp(E)*c*(f-1))-E;
end

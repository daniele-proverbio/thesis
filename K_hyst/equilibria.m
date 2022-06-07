%% Look for stability and phase portrait
% The code below follows the usual steps for linear stability
% An advance would be to use MATCONT
% Change K to see what happens (bistability is destroyed if outside certain
% range); change color map accordingly

%% Author
% Daniele Proverbio, 25/04/2022
% daniele.proverbio@outlook.com
% University of Luxembourg

%% Prepare env
clear; close all; clc;

%% Perform calculation

syms x            % Working with symbolic manipulation (more precise)
K = [0.1,0.192,0.3];          % Basal expression (constant, within accepted range); normally, K = 0.1
c = 0:0.01:4.5;   % Max Production (control parameter)

solutions=[];   % Vector of equilibria solutions
c_vector=[];    % Vector of accepted control parameters
cmap = [];      % Vector of colors


for m = 1:length(c)     % Check all c
    
    f = K(3) + (c(m)*(x)*(x))/(1+(x)*(x)) - (x);     % My equation (vector field)
    soly = vpasolve(f == 0, x);         % Look for roots
    f_prime = diff(f);                  % Estimate derivatives
    
    for n = 1:length(soly)
        if isreal(soly(n))              % Only real roots allowed, obviously
            solutions = [solutions,soly(n)];
            c_vector = [c_vector,c(m)];
            if(vpa(subs(f_prime,x,soly(n))) < 0 )   % Check for linear stability -> eigenvalue <> 0
                cmap = [cmap; [0.5,0.5,1]];             % If Stable point, color it in blue
            else
                cmap = [cmap; [1,0.5,0.5]];             % If ustable, color it in red
            end
        end
    end
end

% plot
figure;
scatter(c_vector,solutions,10,cmap,'filled');
ylim([0 3])
xlim([1, 3])
ax = gca;
ax.FontSize = 18; 
set(gca,'XMinorTick','on','YMinorTick','on')
%title('Equilibria of the system','FontSize',15);
xlabel('c','FontSize',30);
ylabel('$\mathbf{\tilde{x}}$','FontSize',30,'interpreter','latex');
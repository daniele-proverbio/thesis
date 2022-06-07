%% Complement study on focal width
% Analyse data from Dai et al, 2015
% Extract focal width and infer EWS behaviour


clear all
close all
clc

%% Load data
% Data from https://datadryad.org/stash/dataset/doi:10.5061/dryad.k30v3
% (pre-processed)

% NB: swap variables to fit parabola

deterioration = readtable('deterioration.csv'); %Fig A
det_x = deterioration.Var2;
det_y = deterioration.Var1;
cv_det_emp = readtable('cv_det_empirical.csv'); % extracted from figures

diluition  = readtable('diluition.csv'); %Fig B
dil_x = diluition.Var2;
dil_y = diluition.Var1;
cv_dil_emp = readtable('cv_dil_empirical.csv'); % extracted from figures


%% Analysis

fw_det = focal_width(det_x,det_y);
fw_dil = focal_width(dil_x,dil_y);

% Assuming that no correction coefficient is taking place (or assuming that
% we're just making a cmparison and that the coefficients are comparable)
rho_det = 1/fw_det;
rho_dil = 1/fw_dil;

%% Plot corresponding theoretical results
% From x' = p + rho*x^2

%% Fold bifurcation diagram
p = 0:0.01:1;
p1 = 0:0.001:0.1;
x_eq_det_plus = sqrt(p/rho_det);
x_eq_det_minus = - sqrt(p/rho_det);
x_eq_det_plus_s = sqrt(p1/rho_det);
x_eq_det_minus_s = - sqrt(p1/rho_det);

x_eq_dil_plus = sqrt(p/rho_dil);
x_eq_dil_minus = - sqrt(p/rho_dil);

% Conversion to theoretical normal forms
figure()
hold on
plot(p,x_eq_det_plus,color='blue',LineWidth=1)
plot(p,x_eq_det_minus,color='red',LineWidth=1)

plot(p,x_eq_dil_plus,'--',color=[0,154/255,23/255],LineWidth=1)
plot(p,x_eq_dil_minus,'--',color=[1,0.5,0],LineWidth=1)

title('Theoretical normal-form bifurcation diagram')
legend({'Dilution Factor: stable branch','Dilution Factor: unstable branch','Sucrose: stable branch','Sucrose: unstable branch'},Location="northwest",fontsize=13)
ax = gca;
ax.FontSize = 15;
xlabel('p', FontSize=24,Interpreter='latex')
ylabel('$\tilde{x}$', FontSize=24,Interpreter='latex')
hold off

axes('Position',[.6 .2 .3 .3])
box on
hold on
ending = 10;
plot(p1,x_eq_det_plus_s,color='blue',LineWidth=1)
plot(p1,x_eq_det_minus_s,color='red',LineWidth=1)


%% Empirical data
figure(position=[100,100,450,900])

x1 = subplot(2,1,1);
plot(abs(1650 - det_y),det_x,'o',color='blue',LineWidth=1.5)
ax = gca;
ax.FontSize = 16;
ylabel('Population density (cells/$\mu$l)', FontSize=24,Interpreter='latex')
xlabel('DL - Dilution factor', FontSize=24,Interpreter='latex')


x1 = subplot(2,1,2);
plot(abs(0.11 - dil_y),dil_x,'o',color=[0,154/255,23/255],LineWidth=1.5)
ax = gca;
ax.FontSize = 16;
ylabel('Population density (cells/$\mu$l)', FontSize=24,Interpreter='latex')
xlabel('S - [Sucrose] (\%)', FontSize=24,Interpreter='latex')
xlim([0,0.6])

sgtitle('Empirical bifurcation diagrams',fontsize=28,fontweight='bold')

%% CV

sigma = 0.1;  % Random noise intensity (it's unknwn, but anyway it's just a tuning factor)
CV_det = sqrt(sigma*rho_det) ./ (sqrt( p.*sqrt(p) ));
CV_dil = sqrt(sigma*rho_dil) ./ (sqrt( p.*sqrt(p) ));

figure()
hold on
plot(p,CV_det,LineWidth=1.5,color='blue')
plot(p,CV_dil,LineWidth=1.5,color=[0,154/255,23/255])
legend({'Dilution Factor','[Sucrose] (%)'},Location="northeast",fontsize=15)
ylim([0,0.001])
ax = gca;
ax.FontSize = 15;
ylabel('CV', FontSize=18,Interpreter='latex')
xlabel('p', FontSize=18,Interpreter='latex')

sgtitle('Theoretical coefficients of variation',fontsize=18,fontweight='bold')

%%
figure(position=[100,100,450,900])

x1 = subplot(2,1,1);
plot(abs(1650 - cv_det_emp.Var1),cv_det_emp.Var2,'o',color='blue',LineWidth=1.5)
ax = gca;
ax.FontSize = 16;
ylabel('CV', FontSize=24,Interpreter='latex')
xlabel('DL -  Dilution factor', FontSize=24,Interpreter='latex')
ylim([0,0.25])
xlim([0,1000])


x1 = subplot(2,1,2);
plot(abs(0.11 - cv_dil_emp.Var1),cv_dil_emp.Var2,'o',color=[0,154/255,23/255],LineWidth=1.5)
ax = gca;
ax.FontSize = 16;
ylabel('CV', FontSize=24,Interpreter='latex')
xlabel('S - [Sucrose] (\%)', FontSize=24,Interpreter='latex')
xlim([0,0.6])
ylim([0,0.25])

sgtitle('Empirical coefficients of variation',fontsize=28,fontweight='bold')


%% Simulate SDE with a slowly varying parameter
% To see evolution, hysteresis and compute summary statistics measurements
% Model for protein production (from Sharma 2016 or Strogatz
% 2015) as toy model

%% Author
% Daniele Proverbio, 21/01/2022
% daniele.proverbio@outlook.com
% University of Luxembourg


%% NB
% See also "MATLAB/simulate_SDE/normal_form/varying_parameter_NF.m" for
% other plots
%% Prepare env
clear;  clc;
close all;

%% Make Simulation

max_sims = 100;         % how many repeated simulations I want to do
record_tic = zeros(max_sims,1);
record_paramc = zeros(max_sims,1);

% Initialize
time_start = 1;       % Start time of simulations 
time_stop = 1200;     % End recording simulations 
dt = 0.5;             % Time step
sol = zeros((time_stop-time_start),1);
parameters_spanned = zeros(time_stop-time_start,1);

noise = 0.1;           % Noise level (diffusion term) -> pretty low, at least lower than basin height

% SDE simulator
% Euler Maruyama scheme
for jj = 1 : max_sims

    tic = 0;
    f = 0;
    x_in1 = 2.4;            % Initial condition (on upper branch) .  
    x_in2 = 0.1;            % Initial condition (on lower branch) . 
    vals1 = [0.1, 2.7];     % vals1(1) = Basal expression (constant, within accepted range)    % vals1(2) = Max Production (control parameter)  % Moving from right to left 
    vals2 = [0.1, 1.7];     % vals2(1) = Basal expression (constant, within accepted range)    % vals2(2) = Max Production (control parameter)  % Moving from left to right
    for p = time_start:dt:time_stop
        tic = tic + 1;

        if tic == 1
            sol(tic) = x_in1;
            
            parameters_spanned(1) = vals1(2); % to compare evolution of variable vs evolution of parameter
        else
            vals1(2) = vals1(2) - (0.001);   % + 0.001*rand(1) - 0.0005);     % Increase/descrease (+/-) the value of the parameter (depending on forward or backward trend)
            parameters_spanned(tic) = vals1(2);
            f = protein_production(p-dt,sol(tic-1),vals1);
            sol(tic) = sol(tic-1) + f * dt + noise*sqrt(dt)*randn;
        end
        if (record_tic(jj) ==0 && sol(tic) < 1)   % Put 1 to be consistent with theoretical calculation (see §1.1 in report "Normal Forms"). Put 0.5 to be as general as possible (of course results will be slightly delayed)
            record_tic(jj) = ceil(tic*dt);
            record_paramc(jj) = parameters_spanned(tic);
        end
    end

end



%% Plot

% plot time evolution
figure; 
hold on
plot(dt*(time_start-1:30:size(sol(30:end))-1),sol(30:30:end),'o',linewidth=1);  %,'color',[0.9100    0.4100    0.1700]);  % orange color
%plot(dt*(time_start-1:size(sol(30:end))-1),parameters_spanned(30:end));
xlabel('time [arbitrary unit]',fontsize=16);
ylabel('Observed variable',fontsize=16);
%saveas(gcf,'simulation_SDE_1.png');
hold off

% % plot variable vs parameter
% figure; 
% plot(flip(parameters_spanned(30:end)),sol(30:end));
% title('SDE against parameters');
% xlabel('Parameter');
% ylabel('Produced protein');


%% Density function

bins = 25;

%histogram(sol(30:end),bins);
figure
hold on
histogram(sol(30:900),bins,'Normalization','pdf');
histogram(sol(600:1200),bins,'Normalization','pdf');
histogram(sol(1000:1800),bins,'Normalization','pdf');
legend("Time = [0; 900]", "Time = [600; 1200]","Time = [1000; 1800]",fontsize=14)
xlabel('Observed variable values',fontsize=16)
ylabel('Counts (normalised)',fontsize=16)


%% Changepoint analysis

[TF,S1,S2] = ischange(sol(30:30:end)','linear','Threshold',0.3);%,'Threshold',200);
segline = S1.*(1:length(sol(30:30:end))) + S2;
figure
hold on
plot(1:length(sol(30:30:end)),sol(30:30:end),'o',linewidth=1) 
plot(1:length(sol(30:30:end)),segline,linewidth=1)
legend('Data','Linear Regime',fontsize=14,location='northwest')
xlabel("time [arbitrary unit]",fontsize=16)
ylabel("Observed variable",fontsize=16)



%% Equation
% Simulate simple equation for autocatalytic protein production (deterministic part) (Sharma,
% 2015)

function dxdt = protein_production(t,x,vals)

K = vals(1);    
c = vals(2);    

dxdt = x*(1-x/10)-c*(x*x)/(1+x*x); % harvested population   % protein production: K + (c*x*x)/(1+x*x) - x;

end

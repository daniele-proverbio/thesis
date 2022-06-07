%% SDE solver
% Simulate SDE and estimate autocorrelation as summary statistics measurement
% It is connected to the problem of Critical Transitions and Early Warning
% Signals

%% Author
% Daniele Proverbio, 20/08/2019
% daniele.proverbio@outlook.com
% University of Luxembourg


%% Prepare env
clear; close all; clc;

%% Initialize

time_start = 1;     % Start time of simulations 
time_stop = 1000;   % End recording simulations 
dt = 0.5;           % Time step (quite coarse)
sol = zeros(((time_stop-time_start)/dt+1),1);
tic = 0;

x_in = 1.3;             % Initial condition (on upper branch) .  (initial conditions: x_in = [2.25, 1.3])
vals = [0.1, 1.83];     % Basal expression (constant, within accepted range)    % Max Production (control parameter) (values: [2.7, 1.83]) 
noise = 0.05;           % Noise level (diffusion term) -> pretty low, at least lower than basin height

%% SDE simulator
% Euler Maruyama scheme
for p = time_start:dt:time_stop
    tic = tic + 1;
    
    if tic == 1
        sol(tic) = x_in;
    else
        f = protein_production(p-dt,sol(tic-1),vals);
        sol(tic) = sol(tic-1) + f * dt + noise*sqrt(dt)*randn;
    end
    
end

%plot
figure; 
plot(dt*(time_start-1:size(sol)-1),sol(:),'color',[0.9100    0.4100    0.1700]);  % orange color
title('Solution of SDE');
xlabel('time');
ylabel('Produced protein');
%saveas(gcf,'simulation_SDE_1.png');


%% Equation
% Simulate simple equation for autocatalytic protein production (deterministic part) (Sharma,
% 2015)

function dxdt = pitchfork(t,x,vals)

K = vals(1);    
c = vals(2);    

dxdt = K + (c*x*x)/(1+x*x) - x;

end
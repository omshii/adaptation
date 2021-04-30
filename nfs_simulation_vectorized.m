% Plots the average of stochastic simulations over certain time and
% parameters for nfs system. Uses gillespie_vectorized.m (vectorized). 

% TODO: Make into function?

reactants = [0 0]; %A and B

%A production, A degradation, B production, B degradation
reactions = [1 0; -1 0; 0 1; 0 -1]; 

%Parameters
param.I = 0.2; %Input stimulus
param.A0 = 100; 
param.B0 = 100;

%Standard/known parameter set
params = [2 2 10 0.01 4 0.01];

%Assign params
param.k1 = params(1);
param.k2 = params(2);
param.k3 = params(3);
param.K3 = params(4);
param.k4 = params(5);
param.K4 = params(6);

start_time = 0;
end_time = 300; %running time in seconds
dt = 0.1; %in seconds
sims = 200;

%Extend for number of sims
%This is being done here so that you can have different initial values 
reactants = repmat(reactants, sims, 1); 

[time_array1, reactants_array1] = gillespie_vectorized(reactants, reactions, @nfs_propensity_vectorized, param, start_time, 100, dt, sims);
%Increase input stimulus after 1/2 of total time passes
param.I = 0.4;
reactants = reactants_array1(:, :, end);

[time_array2, reactants_array2] = gillespie_vectorized(reactants, reactions, @nfs_propensity_vectorized, param, 100, end_time, dt, sims);

time_array = cat(2, time_array1, time_array2);
reactants_array = cat(3, reactants_array1, reactants_array2);

%Plot A
mean_reactants_array = mean(reactants_array, 1);
mean_reactants_array = reshape(mean_reactants_array, 2, []);
plot(time_array, mean_reactants_array(1, :));
figure(1)

% Random script to test noise_propagation.m
% TODO: clean and delete

%%

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
end_time = 100; %running time in seconds
dt = 0.1; %in seconds
sims = 100;
noise_percent = 25;
reactants = repmat([0 0], sims, 1);

[time_array, complete_traj, input_vector] = noise_propagation(reactants, reactions, @nfs_noisy_propensity_vectorized, param, start_time, end_time, dt, sims, noise_percent);

reactants_traj = complete_traj(:, 1:end-1, :);
input_traj = complete_traj(:, end, :);
%%
input_traj = reshape(input_traj, sims, []);
plot(time_array, input_traj);

%% Test plotting

mr = mean(reactants_array, 1);
ma = reshape(mr, 2, []);
plot(time_array, ma(1, :));

%%
for i=0:4
    figure(i+1)
    ind1 = i*100+1;
    ind2 = (i+1)*100;
    mean_reactants_array = mean(reactants_array(ind1:ind2, :, :), 1);
    mean_reactants_array = reshape(mean_reactants_array, 2, []);
    plot(time_array, mean_reactants_array(1, :));  
end
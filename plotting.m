options=odeset('AbsTol',1e-10,'RelTol',1e-10);

input = [1.05 1 0.95];
noise = [1 5];

I1 = 0.2;

% for i = 1:size(input) 
i = 1;
I2 = 1.05*input(i);
sims = 100;
reactants = [0 0];

start_time = 0;
end_time = 100;
perturb_time = 30;

ffs_params = [10 100 0.1 0.001 1 I1 I2 perturb_time];

[time,proteins] = ode15s(@ffs_ode_mod,[start_time, end_time],reactants,options, ffs_params);
reactants = proteins(end, :);

A_0_index = find(time > 10);
A_0 = proteins(A_0_index(1));

A_1_index = find(time > 90);
A_1 = proteins(A_1_index(1));

delta_A = (abs(A_1 - A_0)/A_0)*100;

delta_I = (abs(I1 - I2)/I1)*100;

susceptibility = delta_A/delta_I;


%nested loop

%for j = 1:size(noise)
j = 1;
reactants = repmat(reactants, sims, 1);

param.k1 = ffs_params(1);
param.k2 = ffs_params(2);
param.k3 = ffs_params(3);
param.K3 = ffs_params(4);
param.k4 = ffs_params(5);
param.I = I2;

noise_percent = noise(j);

[time_array, reactants_array, input_array] = noise_propagation(reactants, reactions, @ffs_propensity_vectorized, param, start_time, perturb_time, dt, sims, noise_percent);

mean_reactants_array = mean(reactants_array, 1);
mean_reactants_array = reshape(mean_reactants_array, 2, []);

mean_A = mean(mean_reactants_array);
std_A = std(mean_reactants_array);

mean_I = mean(input_array);
std_I = std(input_array);

noise_propagation = (std_A/mean_A)/(std_I/mean_I)

options=odeset('AbsTol',1e-10,'RelTol',1e-10);

input = [1.5 1 0.5];
noise = [50 25];
susceptibility = zeros(length(input), 1);
noise_amplification = zeros(length(input), length(noise));

reactants = [0 0]; %A and B

%A production, A degradation, B production, B degradation
reactions = [1 0; -1 0; 0 1; 0 -1]; 


%% Stochastic Simulation

I1 = 0.2;
I2 = 0; %Some default value
perturb_time = 50;
start_time = 0;
end_time = 100;

params = [10 100 0.1 0.001 1 I1 I2 perturb_time];

for i = 1:length(input)
    
    I2 = I1*input(i);
    params = [10 100 0.1 0.001 1 I1 I2 perturb_time];
    
    [time,proteins] = ode15s(@ffs_ode_mod,[start_time, end_time],[0 0],options,params);
    reactants = proteins(end, :);

    A_0_index = find(time < perturb_time);
    A_0 = proteins(A_0_index(end), 1);
    
    if(input(i) < 1)
        peaks = findpeaks(-1*proteins(A_0_index(end):end, 1));
    else
        peaks = findpeaks(proteins(A_0_index(end):end, 1));
    end
    
    A_1 = abs(max(peaks));

    delta_A = (abs(A_1 - A_0)/A_0);
    delta_I = (abs(I1 - I2)/I1);

    susceptibility(i) = delta_A/delta_I;
    
    %figure(i)
    %plot(time, proteins(:, 1))
    
end

%% Deterministic simulation

sims = 100;
dt = 0.1;

param.A0 = 100;
param.B0 = 100;
param.k1 = params(1);
param.k2 = params(2);
param.k3 = params(3);
param.K3 = params(4);
param.k4 = params(5);

start_time = 0;
perturb_time = 100;
end_time = 300;

fig_c = 1;

for i = 1:length(input)
    
    for j = 1:length(noise)
        
        reactants = repmat([0 0], sims, 1);

        noise_percent = noise(j);
        
        param.I = I1*input(i);
        
        [time_array, reactants_array, input_vector] = noise_propagation(reactants, reactions, @ffs_noisy_propensity_vectorized, param, start_time, end_time, dt, sims, noise_percent);
        
        index = find(time_array > perturb_time);
        index = index(1);
        reactants_array = reactants_array(:, :, index:end);
        input_vector = input_vector(index:end);
        
        mean_reactants_array = mean(reactants_array, 1);
        mean_reactants_array = reshape(mean_reactants_array, 2, []);
        
        mean_A = mean(mean_reactants_array(1, :));
        
        std_A = std(mean_reactants_array(1, :));

        mean_I = mean(input_vector);
        std_I = std(input_vector);

        noise_amplification(i, j) = (std_A/mean_A)/(std_I/mean_I);
        
        disp(noise(j));
        disp(input(i));
        disp(noise_amplification(i, j));
    
    end
end

%% Plotting

%scatter(susceptibility', noise_amplification);
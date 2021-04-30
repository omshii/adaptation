% Reads a matrix of nfs parameter values
% Calculates noise amplification and susceptibility
% Writes matrix to file.

options=odeset('AbsTol',1e-10,'RelTol',1e-10);

noise = 25;

reactants = [0 0]; %A and B

%A production, A degradation, B production, B degradation
reactions = [1 0; -1 0; 0 1; 0 -1]; 

I1_array = [0.1 0.2 0.3 0.4];
I2_array = 0.01*I1_array + I1_array;

param_sets = readmatrix('nfs_param_manual.csv');
num_sets = size(param_sets, 1);

susceptibility = zeros(num_sets, length(I1_array));
noise_amplification = zeros(num_sets, length(I1_array));

I1 = 0;
I2 = 0; %Some default value
perturb_time = 50;
start_time = 0;
end_time = 100;

for p=1:num_sets

    for i = 1:length(I1_array)

        I1 = I1_array(i);
        I2 = I2_array(i);

        params = [param_sets(p, :) I1 I2 perturb_time];

        [time,proteins] = ode15s(@nfs_ode_mod,[start_time, end_time],[0 0],options,params);

        reactants = proteins(end, :);

        A_0_index = find(time < perturb_time);
        A_0 = proteins(A_0_index(end), 1);

        peaks = findpeaks(proteins(A_0_index(end):end, 1));

        A_1 = abs(max(peaks));

        delta_A = (abs(A_1 - A_0)/A_0);
        delta_I = (abs(I1 - I2)/I1);

        susceptibility(p, i) = delta_A/delta_I;

    end

    sims = 100;
    dt = 0.1;

    params = param_sets(p, :);
    
    param.A0 = 100;
    param.B0 = 100;
    param.k1 = params(1);
    param.k2 = params(2);
    param.k3 = params(3);
    param.K3 = params(4);
    param.k4 = params(5);
    param.K4 = params(6);

    start_time = 0;
    perturb_time = 20;
    end_time = 100;

    for i = 1:length(I1_array)

            reactants = repmat([0 0], sims, 1);

            noise_percent = noise;

            param.I = I1_array(i);

            [time_array, reactants_array, input_vector] = noise_propagation(reactants, reactions, @nfs_noisy_propensity_vectorized, param, start_time, end_time, dt, sims, noise_percent);

            reactants_array = reactants_array(:, 1:end-1, :);
            
            index = find(time_array > perturb_time);
            index = index(1);
            reactants_array = reactants_array(:, :, index:end);
            input_vector = input_vector(index:end);

            mean_reactants_array = mean(reactants_array, 1);
            mean_reactants_array = reshape(mean_reactants_array, 2, []);

            mean_A = mean(mean_reactants_array(1, :));       
            std_A = std(mean_reactants_array(1, :));

            mean_B = mean(mean_reactants_array(2, :));       
            std_B = std(mean_reactants_array(2, :));

            mean_I = mean(input_vector);
            std_I = std(input_vector);

            noise_amplification(p, i) = (std_A/mean_A)/(std_I/mean_I);
    end
end

plot(susceptibility, noise_amplification, ".", "MarkerSize", 20);
axis([0 1 0 1]);

writematrix(susceptibility, "nfs_susceptibility_1i_25n.csv");
writematrix(noise_amplification, "nfs_noise_1i_25n.csv");
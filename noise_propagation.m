function [time_vector, complete_trajectory, input_vector] = noise_propagation(reactants, reactions, propensity, params, start_time, end_time, dt, sims, noise_percent)
%noise_propagation runs gillespie's algorithm a number of times with noisy input
%
%   reactants: two dimensional array of starting reactant values
%   reactions: array describing reactions
%   propensity: propensity function to be used
%   params: two dimensional array of parameters
%   start_time: simulation starting time
%   end_time: simulation ending time
%   dt: timesteps to use
%   sims: number of simulations
%	noise_percent: percent of noise, for eg. 5 => I between 0.95 I and 1.05 I. 
%
%   time_vector: array of timesteps
%   complete_trajectory: 3D array of sims particle numbers over time

time = repmat(start_time, sims, 1);
time_vector = start_time:dt:end_time;
time_array = repmat(time_vector, sims, 1);

%input_array = zeros(sims, 1, length(time_vector));
%input_trajectory = zeros(sims, 1, length(time_vector));

num_complete = 1;

particle_numbers = zeros(sims, size(reactants, 2), length(time_vector));
complete_trajectory = zeros(sims, size(reactants, 2), length(time_vector));
N = length(time_vector);
counter = ones(sims, 1);

particle_numbers(:, :, 1) = reactants;
counter = counter + 1;

I_lower_bound = 1 - noise_percent/100;
I_upper_bound = 1 + noise_percent/100;

input_vector = ((I_upper_bound-I_lower_bound)*rand(size(time_vector)) + I_lower_bound)*params.I;
input_array = repmat(input_vector, sims, 1);
input = input_vector(1)*ones(sims, 1);

while ~isempty(time)
    
    v = propensity(reactants, params, input);
    h = sum(v, 2);
    
    %Generate 2 random numbers
    r1 = rand(sims, 1);
    r2 = rand(sims, 1);

    tau = -log(r1)./h; %calculate time interval
    nv = bsxfun(@rdivide, v, h); %normalize reaction rates
    cl = cumsum(nv, 2); %cumulative sum of rates
    
    %find reaction to occur next
    index = length(reactions) - sum(bsxfun(@gt, cl, r2), 2) + 1;

    %update reactants
    reactants = reactants + reactions(index, :);
    
    %update time
    time = time + tau;
    
    update = time > time_array(sub2ind(size(time_array), 1:sims, counter'))';
    
    while any(update)
        particle_numbers(sub2ind(size(particle_numbers), 1:sims, ones(1, sims), counter')) = reactants(:,1) .* update;
        particle_numbers(sub2ind(size(particle_numbers), 1:sims, ones(1, sims)*2, counter')) = reactants(:,2) .* update;
        input = ones(size(input)).*(input_array(counter));
        counter = counter+update;
        completed = counter > N;
        if any(completed)
            complete_trajectory(num_complete:num_complete+sum(completed)-1, :, :) = particle_numbers(completed, :, :);
            num_complete = num_complete+sum(completed);
            counter(completed) = [];
            time(completed) = [];
            time_array(completed, :) = [];
            reactants(completed, :) = [];
            input(completed) = [];
            input_array(completed, :) = [];
            particle_numbers(completed, :, :) = [];
            sims = sims - sum(completed);
            if isempty(counter)
                break
            end
        end
        update = time > time_array(sub2ind(size(time_array), 1:sims, counter'))';
    end
    
end

end

%TODO: Add vectorization for different parameter values


function [time_vector, complete_trajectory] = gillespie_vectorized(reactants, reactions, propensity, params, start_time, end_time, dt, sims)

reactants = repmat(reactants, sims, 1);
time = repmat(start_time, sims, 1);
time_vector = start_time:dt:end_time;
time_array = repmat(time_vector, sims, 1);
num_complete = 1;

particle_numbers = zeros(sims, size(reactants, 2), length(time_vector));
complete_trajectory = zeros(sims, size(reactants, 2), length(time_vector));
N = length(time_vector);
counter = ones(sims, 1);

particle_numbers(:, :, 1) = reactants;
counter = counter + 1;

while ~isempty(time)
    
    v = propensity(reactants, params);
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
        counter = counter+update;
        completed = counter > N;
        if any(completed)
            complete_trajectory(num_complete:num_complete+sum(completed)-1, :, :) = particle_numbers(completed, :, :);
            num_complete = num_complete+sum(completed);
            counter(completed) = [];
            time(completed) = [];
            time_array(completed, :) = [];
            reactants(completed, :) = [];
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


function [time_array, particle_numbers] = gillespie(reactants, reactions, propensity, params, start_time, end_time, dt)

time = start_time;
time_array = start_time:dt:end_time;
particle_numbers = zeros(length(time_array), length(reactants));
N = length(time_array);
counter = 1;

particle_numbers(1, :) = reactants;
counter = counter + 1;

while time < time_array(end) && counter <= N
    
    v = propensity(reactants, params);
    h = sum(v);
    
    %Generate 2 random numbers
    r1 = rand();
    r2 = rand();

    tau = -log(r1)/h; %calculate time interval
    nv = v/h; %normalize reaction rates
    cl = cumsum(nv); %cumulative sum of rates
    
    %find reaction to occur next
    temp = find(cl > r2);
    index = temp(1);

    %update reactants
    reactants = reactants + reactions(index, :);
    
    %update time
    time = time + tau;
    
    while counter<=N && time>time_array(counter)
        particle_numbers(counter, :) = reactants(:);
        counter = counter+1;
    end
    
end

end


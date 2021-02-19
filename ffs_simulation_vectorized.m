reactants = [0 0]; %A and B

%A production, A degradation, B production, B degradation
reactions = [1 0; -1 0; 0 1; 0 -1]; 

%Parameters
param.I = 0.2; %Input stimulus
param.A0 = 100; 
param.B0 = 100;

%Standard/known parameter set
params = [10 100 0.1 0.001 1];

%Filtered params
params = [1.01175338050600,28.4905059445019,0.109854722130800,0.846454413958243,0.101186727319175]; %TRYME
%params = [0.503208702911553,87.6963600687840,8.30251100085423,0.00158775170497156,0.821904001529459];
%params = [45.0391467316844,637.044241305400,5.88014725535334,0.00295193612580742,0.167903768689519];

%Assign params
param.k1 = params(1);
param.k2 = params(2);
param.k3 = params(3);
param.K3 = params(4);
param.k4 = params(5);

start_time = 0;
end_time = 1000; %running time in seconds
dt = 0.1; %in seconds
sims = 500;
reactants = repmat(reactants, sims, 1);

[time_array1, reactants_array1] = gillespie_vectorized(reactants, reactions, @ffs_propensity_vectorized, param, start_time, 500, dt, sims);
%Increase input stimulus after 1/2 of total time passes
param.I = 0.4;
reactants = reactants_array1(:, :, end);

[time_array2, reactants_array2] = gillespie_vectorized(reactants, reactions, @ffs_propensity_vectorized, param, 500, end_time, dt, sims);

time_array = cat(2, time_array1, time_array2);
reactants_array = cat(3, reactants_array1, reactants_array2);

%Plot A
mean_reactants_array = mean(reactants_array, 1);
mean_reactants_array = reshape(mean_reactants_array, 2, []);
plot(time_array, mean_reactants_array(1, :));
figure(1)

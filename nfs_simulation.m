%Simulation function, single run, non vectorized

reactants = [0 0]; %A and B

%A production, A degradation, B production, B degradation
reactions = [1 0; -1 0; 0 1; 0 -1]; 

%Parameters
param.I = 0.2; %Input stimulus
param.A0 = 100; 
param.B0 = 100; 
param.k1 = 2;
param.k2 = 2;
param.k3 = 10;
param.K3 = 0.01;
param.k4 = 4;
param.K4 = 0.01;

start_time = 0;
end_time = 50; %running time in seconds
dt = 0.1; %in seconds

[X, Y] = gillespie(reactants, reactions, @nfs_propensity, param, start_time, end_time, dt);

figure(1)
plot(X, Y(:, 1))

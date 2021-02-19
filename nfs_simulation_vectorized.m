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

params = [58.9002521426813,0.226012097687259,1.24316313998002,0.00388926446305617,0.394658563388276,0.00453088748125769];
params = [21.6739748798593,0.267154330430257,15.8254716453190,0.00105520460572554,6.07201675776149,0.0301616987181525];
params = [35.6493871070062,2.13078712696238,6.52528986565635,0.00120582618146524,2.48245486491569,0.00106420209362253];
params = [7.56256982996279,0.552844674236217,0.629574034686460,0.00155861231247859,0.241346919433197,0.00648982588403310];
%params = [47.5082469452525,0.144106817638821,1.62641666825525,0.00220252001007801,0.335651786238324,0.143680144503080];
%params = [0.105999877847376,5.13087094530508,12.0663336616234,0.00470438477559106,13.8245685509548,0.0526934345787161];
%params = [16.1001, 2.9939, 0.2810, 0.0609, 44.1849, 0.0096, 1.5439];
%params = [0.272332581610224,1.50368681241021,7.21352692794849,0.0255900617560167,10.0619553285569,0.329216385829335];
%params = [3.48673558163433,2.06170320236162,77.9466373920635,0.00464975668282479,6.07618955873295,0.0236613606258170];
%params = [8.16284358024069,40.2663471990195,50.5476976378181,0.217404776960704,23.8375595794157,0.0382645315601751];
%params = [8.99312641566308,5.01358782249180,52.5363745816817,0.102651412192301,0.661190857352090,0.0198810339184314];
params =[0.828260904155780,24.5471265048571,14.8265959483464,0.0194773931397289,0.244079592636333,0.00451453608258908] ;
params = [0.828260904155780,24.5471265048571,14.8265959483464,0.0194773931397289,0.244079592636333,0.00451453608258908];

param.k1 = params(1);
param.k2 = params(2);
param.k3 = params(3);
param.K3 = params(4);
param.k4 = params(5);
param.K4 = params(6);

start_time = 0;
end_time = 50; %running time in seconds
dt = 0.1; %in seconds
sims = 500;
reactants = repmat(reactants, sims, 1);

[time_array1, reactants_array1] = gillespie_vectorized(reactants, reactions, @nfs_propensity_vectorized, param, start_time, 14, dt, sims);
%Increase input stimulus after 1/4th of total time passes
param.I = 0.4;
reactants = reactants_array1(:, :, end);

[time_array2, reactants_array2] = gillespie_vectorized(reactants, reactions, @nfs_propensity_vectorized, param, 14, end_time, dt, sims);

time_array = cat(2, time_array1, time_array2);
reactants_array = cat(3, reactants_array1, reactants_array2);

%figure(1)
mean_reactants_array = mean(reactants_array, 1);
mean_reactants_array = reshape(mean_reactants_array, 2, []);
plot(time_array, mean_reactants_array(1, :));
figure(1)

%Bm = mean(B, 1);
%Bc = reshape(Bm, 2, []);
%plot(A, Bc(1, :));
%figure(2)y
%plot(A, Bc(2, :));

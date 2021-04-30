% Script to submit random_parameters.m 
% to Colgate's Turing cluster
% Generate and store random parameter sets 

c1 = parcluster('Turing');
j1 = batch(c1, @random_parameters, 4, {}, 'Pool', 3);
wait(j1);
output = fetchOutputs(j1);
save("output.mat", "output");
csvwrite("nfs_params", output{1, 2});
csvwrite("ffs_params", output{1, 4});
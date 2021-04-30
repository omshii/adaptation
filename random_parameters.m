function [nfs_num_params, nfs_params, ffs_num_params, ffs_params] = random_parameters()
%random_parameters: generates X random parameters for nfs and ffs systems
%                   ,filters them using filter_parameters, and writes to
%                   file.

%TODO: Add X as parameter to function
%      Add X to filenames

X = 10000;
str_stamp = erase(string(now), ".");

%For NFS
k1 = -1 + (2+1)*rand(X,1);
k2 = -1 + (2+1)*rand(X,1);
k3 = -1 + (2+1)*rand(X,1);
k4 = -1 + (2+1)*rand(X,1);
K3 = -3 + (0+3)*rand(X,1);
K4 = -3 + (0+3)*rand(X,1);

k1 = 10.^k1;
k2 = 10.^k2;
k3 = 10.^k3;
k4 = 10.^k4;
K3 = 10.^K3;
K4 = 10.^K4;

nfs_params = [k1 k2 k3 K3 k4 K4];
N = size(nfs_params, 2);

[nfs_num_params, nfs_params] = filter_params(nfs_params, @nfs_ode);
nfs_filtered_params = nfs_params((nfs_params(:, N+1)>=1 & nfs_params(:, N+2)>=10), :);
save(str_stamp+"_nfs_params.mat", "nfs_params");
save(str_stamp+"_nfs_params_filtered.mat", "nfs_filtered_params");

disp(nfs_num_params);



%For FFS
k1 = -1 + (2+1)*rand(X,1);
k2 = -1 + (3+1)*rand(X,1);
k3 = -1 + (2+1)*rand(X,1);
k4 = -1 + (2+1)*rand(X,1);
K3 = -3 + (0+3)*rand(X,1);

k1 = 10.^k1;
k2 = 10.^k2;
k3 = 10.^k3;
k4 = 10.^k4;
K3 = 10.^K3;

ffs_params = [k1 k2 k3 K3 k4];
N = size(ffs_params, 2);

[ffs_num_params, ffs_params] = filter_params(ffs_params, @ffs_ode);
ffs_filtered_params = ffs_params((ffs_params(:, N+1)>=1 & ffs_params(:, N+2)>=10), :);
save(str_stamp+"_ffs_params.mat", "ffs_params");
save(str_stamp+"_ffs_params_filtered.mat", "ffs_filtered_params");

disp(ffs_num_params);

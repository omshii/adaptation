function [nfs_num_params, nfs_params, ffs_num_params, ffs_params] = random_parameters()

X = 10000;

% Generate X number of random parameter sets, and call filter_params.

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

[nfs_num_params, nfs_params] = filter_params(nfs_params, @nfs_ode);

%disp(nfs_num_params);

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

[ffs_num_params, ffs_params] = filter_params(ffs_params, @ffs_ode);

%disp(ffs_num_params);





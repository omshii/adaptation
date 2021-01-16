X = 100;

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

params = [k1 k2 k3 K3 k4 K4];

num_params = filter_params(params, @nfs_ode);

disp(num_params);


%TODO: Add ffs param filtering.




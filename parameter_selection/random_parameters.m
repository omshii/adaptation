% Generate 2 100 random parameter sets and save to file

%For NFS
k1 = -1 + (2+1)*rand(100,1);
k2 = -1 + (2+1)*rand(100,1);
k3 = -1 + (2+1)*rand(100,1);
k4 = -1 + (2+1)*rand(100,1);
K3 = -3 + (0+3)*rand(100,1);
K4 = -3 + (0+3)*rand(100,1);

k1 = 10.^k1;
k2 = 10.^k2;
k3 = 10.^k3;
k4 = 10.^k4;
K3 = 10.^K3;
K4 = 10.^K4;

params = [k1 k2 k3 K3 k4 K4];






output = matfile("output.mat");
nfs_params = output.output(1, 2);
ffs_params = output.output(1, 4);

nfs_params = nfs_params{1, 1};
ffs_params = ffs_params{1, 1};

ffs_filtered_params = ffs_params((ffs_params(:, 7)>=1 & ffs_params(:, 8)>=10), :);
nfs_filtered_params = nfs_params((nfs_params(:, 6)>=1 & nfs_params(:, 7)>=10), :);


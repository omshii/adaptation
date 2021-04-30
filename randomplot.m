% Reading ffs/nfs noise propagation vs susceptibility files and plotting
% TODO: Move this to ffs/nfs_plotting.m

nfs_susceptibility = reshape(readmatrix("nfs_susceptibility_1i_10n.csv"), [], 1);
nfs_noise = reshape(readmatrix("nfs_noise_1i_10n.csv"), [], 1);

plot(nfs_susceptibility, nfs_noise, '.k', 'MarkerSize', 20);

hold on

ffs_susceptibility = reshape(readmatrix("ffs_susceptibility_1i_10n.csv"), [], 1);
ffs_noise = reshape(readmatrix("ffs_noise_1i_10n.csv"), [], 1);

plot(ffs_susceptibility, ffs_noise, '.r', 'MarkerSize', 20);
hold on

legend('negative feedback', 'feed forward');

xlabel("susceptibility");
ylabel("noise amplification");
axis([0 1.5 0 1.5])
print("ffs_nfs_1i_10n", "-dpng");
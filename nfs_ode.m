function dydt = nfs_ode(t, y, params)
%nfs_ode describes odes for nfs
%
%   t: the current time
%   y: set of current reactant values, of the form [A B]
%   params: set of parameters containing-
%            k1, k2, k3, K3, k4, K4
%
%   dydt: set of calculated derivatives of the form [dA/dt dB/dt]

A = y(1);
B = y(2);

if(t > 500)
    I = 0.4;
else
    I = 0.2;
end

dA = params(1)*I*(1-A)-params(2)*A*B;
dB = params(3)*A*(1-B)/(params(4)+1-B)-params(5)*B/(params(6)+B);

dydt = [dA; dB];
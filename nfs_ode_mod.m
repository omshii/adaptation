function dydt = nfs_ode_mod(t, y, params)
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

I = params(7);

if(t > params(9))
    I = params(8);
end

dA = params(1)*I*(1-A)-params(2)*A*B;
dB = params(3)*A*(1-B)/(params(4)+1-B)-params(5)*B/(params(6)+B);

% pr.k1 = params(1);
% pr.k2 = params(2);
% pr.k3 = params(3);
% pr.K3 = params(4);
% pr.k4 = params(5);
% pr.K4 = params(6);
% pr.A0 = 100;
% pr.B0 = 100;
% pr.I = 0.2;
% 
% dA=pr.k1*pr.I*(pr.A0-A) - (pr.k2*A*B)/pr.B0; 
% dB=pr.k3*A*(pr.B0-B)/(pr.A0*(pr.K3+1-(B/pr.B0))) - (pr.k4*B*pr.B0)/(pr.K4*pr.B0+B); %B degradation 


dydt = [dA; dB];
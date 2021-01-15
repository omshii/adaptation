%NFS Model

function dydt = nfs_ode(t, y, k1, k2, k3, K3, k4, K4)

A = y(1);
B = y(2);

if(t > 10) 
    I = 0.4;
else
    I = 0.2;
end

dA = k1*I*(1-A)-k2*A*B;
dB = k3*A*(1-B)/(K3+1-B)-k4*B/(K4+B);
dydt = [dA; dB];
%NFS Model

function dydt = nfs_ode(t, y, params)

A = y(1);
B = y(2);

if(t > 10) 
    I = 0.4;
else
    I = 0.2;
end

dA = params(1)*I*(1-A)-params(2)*A*B;
dB = params(3)*A*(1-B)/(params(4)+1-B)-params(5)*B/(params(6)+B);
dydt = [dA; dB];
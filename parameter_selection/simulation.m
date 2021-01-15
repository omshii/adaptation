%simulation to filter parameters

k1 = -1 + (2+1)*rand();
k2 = -1 + (2+1)*rand();
k3 = -1 + (2+1)*rand();
k4 = -1 + (2+1)*rand();
K3 = -3 + (0+3)*rand();
K4 = -3 + (0+3)*rand();

k1 = 10^k1;
k2 = 10^k2;
k3 = 10^k3;
k4 = 10^k4;
K3 = 10^K3;
K4 = 10^K4;

%Standard/known params
k1= 2;
k2= 2;
k3= 10;
K3= 0.01;
k4= 4;
K4= 0.01;

start_time = 0;
end_time = 500;

A_0= 1*0.3931;
B_0= 1*0.3088;

[time,proteins]=ode15s(@nfs_ode,[start_time, end_time],[A_0 B_0],[],k1,k2,k3,K3,k4,K4);

slope = (proteins(end, 1) - proteins(end-1, 1))/(time(end, 1) - time(end-1, 1));

peaks = findpeaks(proteins(:, 1));
peaks = abs(peaks - A_0);

num_peaks = length(peaks);
[max_peak, max_index] = max(peaks);

if slope < 0.0001
    %Calculate regular sensitivity and precision
    O_peak = max_peak;
    sensitivity = abs(((O_peak - A_0)/A_0)/(1));
    precision = abs((1)/(abs(A_0 - proteins(end, 1))/A_0));
elseif (num_peaks > 1) && (max_index ~= length(peaks))
    O_peak1 = max_peak;
    O_peak2 = peaks(max_index + 1);
    if(O_peak1 > 2*O_peak2)
        %Calculate oscillating sensitivity and precision
        sensitivity = abs(((O_peak - A_0)/A_0)/(1));
        precision = abs((1)/(abs(A_0 - proteins(end, 1))/A_0));
    else
        pass = 0;
        return
    end
else
    pass = 0;
    return
end 


%Plot for debugging

%figure(1);
%plot(time, proteins(:, 1));
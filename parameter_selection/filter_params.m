function [num_params] = filter_params(params, ode_func)
%filter_params: Checks sensitivity and precision for parameter sets. 
%   ode_func - function describing ode's for each set of parameters.
%   Takes in M X N array - params, which is the random parameter sets. 
%   Returns num_params - number of params found to be adapting.
%   Writes params, with two additional columns for precision and sensitivity to file.  

M = size(params, 1);
N = size(params, 2);
sensitivity = zeros(M, 1);
precision = zeros(M, 1);

A_0= 1*0.3931;
B_0= 1*0.3088;
start_time = 0;
end_time = 500;

for i = 1:M
    
    [time,proteins]=ode15s(ode_func,[start_time, end_time],[A_0 B_0],[],params(i, :));
    
    slope = (proteins(end, 1) - proteins(end-1, 1))/(time(end, 1) - time(end-1, 1));
    
    peaks = findpeaks(proteins(:, 1));
    peaks = abs(peaks - A_0);
    
    num_peaks = length(peaks);
    [max_peak, max_index] = max(peaks);
    
    if slope < 0.0001
        
        %Calculate regular sensitivity and precision
        O_peak = max_peak;
        
        sensitivity(i, 1) = abs(((O_peak - A_0)/A_0)/(1));
        precision(i, 1) = abs((1)/(abs(A_0 - proteins(end, 1))/A_0));
        
    elseif (num_peaks > 1) && (max_index ~= length(peaks))
        
        O_peak1 = max_peak;
        O_peak2 = peaks(max_index + 1);
        
        if(O_peak1 > 2*O_peak2)
            
            %Calculate oscillating sensitivity and precision
            sensitivity(i, 1) = abs(((O_peak - A_0)/A_0)/(1));
            precision(i, 1) = abs((1)/(abs(A_0 - proteins(end, 1))/A_0));
            
        else
            
            sensitivity(i, 1) = -1;
            precision(i, 1) = -1;
            
        end
        
    else
        
        sensitivity(i, 1) = -1;
        precision(i, 1) = -1;
        
    end
end

params(:, N+1) = sensitivity;
params(:, N+2) = precision;

num_params = length(params(params(:, N+1)>1 & params(:, N+2)>10));








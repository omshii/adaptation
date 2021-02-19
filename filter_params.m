function [num_params, params] = filter_params(params, ode_func)
%filter_params checks adaptation sensitivity and precision for parameter sets.
%
%   params: array of N parameter sets, with M elements in each set.
%   ode_func: describes odes on which params to be tested.
%
%   num_params: number of parameters in params satisfying
%               sensitivity > 1 and precision > 10
%   params: params array with added columns for
%               sensitivity, precision, damped

M = size(params, 1);
N = size(params, 2);
sensitivity = zeros(M, 1);
precision = zeros(M, 1);

% is one if oscillations are damped
% (satisfying Opeak1 > 2*Opeak2)
% else zero
damped = zeros(M, 1); 

I_1= 0.2;
I_2= 0.4;
start_time = 0;
end_time = 1000;

parfor i = 1:M

    [time,proteins] = ode15s(ode_func,[start_time, end_time],[0 0],[],params(i, :));
    
    A_0_index = find(time > 450);
    A_0 = proteins(A_0_index(1));
    
    %Indices to check slopes at
    index1 = find(time > 800);
    index1 = index1(1);
    
    index2 = find(time > 900);
    index2 = index2(1);
    
    index3 = length(time);
    
    %Check slope
    slope1 = (proteins(index1, 1) - proteins(index1-1, 1))/(time(index1, 1) - time(index1-1, 1));
    slope2 = (proteins(index2, 1) - proteins(index2-1, 1))/(time(index2, 1) - time(index2-1, 1));
    slope3 = (proteins(index3, 1) - proteins(index3-1, 1))/(time(index3, 1) - time(index3-1, 1));
    
    %Find peaks
    peaks = findpeaks(proteins(:, 1));
    num_peaks = length(peaks); 
    
    if all([slope1 slope2 slope3] < 0.0000001)
     
    %if true
        
        %Calculate sensitivity and precision
        O_peak = max(peaks);
        
        sensitivity(i, 1) = abs(((O_peak - A_0)/A_0)/((I_2 - I_1)/I_1));
        precision(i, 1) = abs((1)/(abs(A_0 - proteins(end, 1))/A_0));
        
    else        
        sensitivity(i, 1) = -1;
        precision(i, 1) = -1;
    end   
    
    if (num_peaks > 1)
        
        [~, trough_indices] = findpeaks(-1*proteins(:, 1));
        [max_peak, Opeak1_index] = max(peaks);
        Opeak2_index = find(trough_indices > Opeak1_index);
        Opeak2_index = Opeak2_index(1);
        
        O_peak1 = abs(A_0 - max_peak);
        O_peak2 = abs(A_0 - proteins(Opeak2_index, 1)); 
        
        if(O_peak1 > 2*O_peak2)
            damped(i, 1) = 1;
        end
    end
            
    
end

params(:, N+1) = sensitivity;
params(:, N+2) = precision;
params(:, N+3) = damped;

num_params = length(params(params(:, N+1)>=1 & params(:, N+2)>=10));

function [x] = eval_movpeaks(gen)
global current_evals,global global_max,global current_error, global recent_change,global change_frequency,global number_of_peaks, global evals,  
if ((change_frequency > 0)&&(mod(evals,change_frequency) == 0))
    change_peaks;
end
dummy = zeros(1,number_of_peaks);
for i = 1:number_of_peaks
    dummy(i) = Peak_Function_Cone(gen , i);
end
evals=evals+1;
current_evals = current_evals + 1;
x = max(dummy);
if (recent_change || (x > (global_max-current_error)))
    current_error = global_max - x;
    recent_change = false;
end
end
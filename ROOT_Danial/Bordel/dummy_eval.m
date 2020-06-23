function [x] = dummy_eval(gen)
global number_of_peaks;
dummy = zeros(1,number_of_peaks);
for i = 1:number_of_peaks
    dummy(i) = Peak_Function_Cone(gen , i);
end
x = max(dummy);
end
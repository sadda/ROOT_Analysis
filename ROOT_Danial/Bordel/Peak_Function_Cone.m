function [x] = Peak_Function_Cone (gen,peak_number) 
global geno_size, global peak; 
    x = peak(peak_number,geno_size+2)-((peak(peak_number,geno_size+1) * sqrt(sum ((gen(1:geno_size)-peak(peak_number,1:geno_size)).^2))));
end
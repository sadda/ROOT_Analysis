function [] = change_peaks()
global current_evals, global recent_change_PSO, global recent_change, global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global minwidth,global maxwidth,  global maximum_peak,  global global_max,global peak,global prev_movement,
for i = 1:number_of_peaks
    R = rand(1,geno_size)-0.5;
    L=1-lambda;
    Shift = vlength(i)  * (((L * R) + (lambda * prev_movement(i,:))) / Euc_dis(((L * R) + (lambda * prev_movement(i,:))) ,zeros(1,geno_size)));
    for j = 1:geno_size
        if ((peak(i,j) + Shift(j)) < mincoordinate)
            dummy =  (2*mincoordinate) - peak(i,j) - Shift(j);
            prev_movement(i,j) =  dummy - peak(i,j);
            peak(i,j) = dummy;
        elseif ((peak(i,j) + Shift(j)) > maxcoordinate)
            dummy =	(2*maxcoordinate) - peak(i,j)	- Shift(j);
            prev_movement(i,j) = dummy - peak(i,j);
            peak(i,j) = dummy;
        else
            peak(i,j) = peak(i,j) + Shift(j);
            prev_movement(i,j) = Shift(j);
        end
    end
    
    %/* change peak width */
    j = geno_size+1;
    offset = randn * width_severity(i);
    if ((peak(i,j) + offset) < minwidth)
        peak(i,j) = 2.0 * minwidth - peak(i,j) - offset;
    elseif ((peak(i,j) + offset) > maxwidth)
        peak(i,j) = 2.0 * maxwidth - peak(i,j) - offset;
    else
        peak(i,j) = peak(i,j) + offset;
    end
    %/* change peak height */
    j=j+1;
    offset = height_severity(i) * randn;
    if ((peak(i,j) + offset) < minheight)
        peak(i,j) = 2.0 * minheight - peak(i,j) - offset;
    elseif ((peak(i,j) + offset) > maxheight)
        peak(i,j) = 2.0 * maxheight - peak(i,j) - offset;
    else
        peak(i,j) = peak(i,j) + offset;
    end
end
%obtaining global max
global_max = -inf;
for i = 1: number_of_peaks
    dummy = dummy_eval(peak(i,1:geno_size));
    if (dummy > global_max)
        global_max = dummy;
        maximum_peak = i;
    end
end
recent_change = true;
recent_change_PSO = 1;
current_evals = 0;
end

function [] =init_peaks()
global maxcoordinate,global mincoordinate,global geno_size,global number_of_peaks,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,   global global_max,global peak,global shift,global prev_movement,
shift = zeros(geno_size,1);
peak = zeros(number_of_peaks ,geno_size + 2);
peak(:,1:geno_size) = mincoordinate + ((maxcoordinate-mincoordinate) * rand(number_of_peaks ,geno_size));
prev_movement = rand(number_of_peaks ,geno_size) - 0.5;

if (standardwidth == 0.0)
    peak(1:number_of_peaks,geno_size+1) = minwidth + ((maxwidth - minwidth) * rand(1,number_of_peaks));
else
    peak(1:number_of_peaks,geno_size+1) = standardwidth * ones(number_of_peaks,1);
end

if (standardheight <= 0.0)
    peak(1:number_of_peaks,geno_size + 2) = minheight + ((maxheight - minheight) * rand(number_of_peaks,1));
else
    peak(1:number_of_peaks,geno_size + 2) = standardheight * ones(number_of_peaks,1);
end
end

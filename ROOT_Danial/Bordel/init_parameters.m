    % number of evaluations between changes. change_frequency
    % =0 means that function never changes (or only if function change_peaks is called)
    % Scenario 1: 5000
    global current_evals, global Delta_e, global recent_change_PSO, global environment_change, global standardwidth;global minwidth;global maxwidth;global minheight;global maxheight;global mincoordinate;global maxcoordinate;global calculate_right_peak;global calculate_offline_performance;global calculate_average_error;global use_basis_function;global lambda;global width_severity;global height_severity;global vlength;global change_frequency;global number_of_peaks;global movrandseed;global geno_size;global standardheight;
    change_frequency = 2500;
    number_of_peaks = 20; % number of peaks in the landscape
    geno_size = 5; % number of dimensions, or the number of  valued genes
    vlength = 0.5 + (2.5 * rand(number_of_peaks,1)); %Shift severity
%     vlength = ones(number_of_peaks,1); %Shift severity
%     threshold = 50;%survival time threshold
    mincoordinate = -50.0;% minimum coordinate in each dimension
    maxcoordinate = 50.0;% maximum coordinate in each dimension
    movrandseed = 1; % seed for built-in random number generator
    height_severity = (14*rand(number_of_peaks,1))+1;% severity of height changes, larger numbers mean larger severity==> range[1-15]
    width_severity = (1.4*rand(number_of_peaks,1))+0.1;% severity of width changes, larger numbers mean larger severity==> range[0.1,1.5]
    lambda = 0;%For lambda = 1.0 each move has the same direction,while for lambda = 0.0, each move has a random direction
    calculate_right_peak = true; 
    minheight = 30.0;% minimum height of the peaks
    maxheight = 70.0;% maximum height of the peaks
    standardheight = 50;% height chosen randomly when standardheight = 0.0
    minwidth = 1;% minimum width of the peaks
    maxwidth =12;% maximum width of the peaks
    standardwidth = 6;% width chosen randomly when standardwidth = 0.0
    environment_change = 100;
    Delta_e = change_frequency-1;
                        %***%END OF PARAMETER SECTION ****%
    global current_error;global evals;
    current_error = inf; % error of the currently best individual %
    evals = 0;
    current_evals = 0;
    recent_change_PSO = 0;
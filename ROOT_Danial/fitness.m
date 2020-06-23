function Pso = fitness(x_in , Pso)
global Delta_e,global current_evals, global threshold, global ROOT;global Last,global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global slope_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minslope,global maxslope, global standardslope,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType; %#ok<NUSED>
[l,~]=size(x_in);
Pso.result=zeros(l,1);
if evals<Pso.max_fe
    for i=1:l
        x=x_in(i,:);
        Pso.result(i,1)=0;
        Pso.result(i,1) = eval_movpeaks(x);
        if ceil((evals+1)/change_frequency)==1
            Last.error(ceil((evals+1)/change_frequency) , Pso.test) = current_error;
            Last.current_max(ceil((evals+1)/change_frequency) , Pso.test) = global_max - current_error;
            Last.global_max(ceil((evals+1)/change_frequency) , Pso.test) = global_max;
        end
        if (current_evals==Delta_e)
            [~, Best_found_index]= max(Pso.Gbest_value);
            %% Choosing the first RS
            if evals<change_frequency%RS for the 1st environment
                ROOT.solution_counter_T40 = 1;
                ROOT.solution_position_T40 = Pso.Gbest_position(Best_found_index,:);
                ROOT.solution_value_T40 = Pso.Gbest_value(Best_found_index);
                ROOT.survival_time_sequence_T40  = 0;
%__________________________________________________________________________________                
                ROOT.solution_counter_T45  = 1;
                ROOT.solution_position_T45  = Pso.Gbest_position(Best_found_index,:);
                ROOT.solution_value_T45  = Pso.Gbest_value(Best_found_index);
                ROOT.survival_time_sequence_T45  = 0;
%__________________________________________________________________________________                
                ROOT.solution_counter_T50  = 1;
                ROOT.solution_position_T50  = Pso.Gbest_position(Best_found_index,:);
                ROOT.solution_value_T50  = Pso.Gbest_value(Best_found_index);
                ROOT.survival_time_sequence_T50  = 0;
%__________________________________________________________________________________                
            else                
                %%  T40
                if ROOT.solution_value_T40  < 40
                    clear dummy;
                    clear RS_index;
                    if sum(40<Pso.Gbest_value-Pso.Fit_var)~=0
                        dummy = 40 + Pso.Fit_var;
                        dummy(find(dummy>Pso.Gbest_value))=inf;
                        dummy(find(dummy<inf))=0;
                        [~,RS_index] = min (dummy + (Pso.vlength./max(Pso.vlength))+(Pso.heigth_var./max(Pso.heigth_var)));
                    else
                        [~,RS_index] = max(Pso.Gbest_value);
                    end
                    ROOT.solution_counter_T40  = ROOT.solution_counter_T40  + 1;
                    ROOT.solution_position_T40  = Pso.Gbest_position(RS_index,:);
                    ROOT.solution_value_T40  = Pso.Gbest_value(RS_index);
                    ROOT.survival_time_sequence_T40 (1,ceil((evals)/change_frequency))=0;
                else
                    ROOT.survival_time_sequence_T40 (1,ceil((evals)/change_frequency)) = ROOT.survival_time_sequence_T40 (1,((ceil((evals)/change_frequency))-1)) + 1;
                end
                %%  T45
                if ROOT.solution_value_T45  < 45
                    clear dummy;
                    clear RS_index;
                    if sum(45<Pso.Gbest_value-Pso.Fit_var)~=0
                        dummy = 45 + Pso.Fit_var;
                        dummy(find(dummy>Pso.Gbest_value))=inf;
                        dummy(find(dummy<inf))=0;
                        [~,RS_index] = min (dummy + (Pso.vlength./max(Pso.vlength))+(Pso.heigth_var./max(Pso.heigth_var)));
                    else
                        [~,RS_index] = max(Pso.Gbest_value);
                    end
                    ROOT.solution_counter_T45  = ROOT.solution_counter_T45  + 1;
                    ROOT.solution_position_T45  = Pso.Gbest_position(RS_index,:);
                    ROOT.solution_value_T45  = Pso.Gbest_value(RS_index);
                    ROOT.survival_time_sequence_T45 (1,ceil((evals)/change_frequency))=0;
                else
                    ROOT.survival_time_sequence_T45 (1,ceil((evals)/change_frequency)) = ROOT.survival_time_sequence_T45 (1,((ceil((evals)/change_frequency))-1)) + 1;
                end
                %% T50
                if ROOT.solution_value_T50  < 50
                    clear dummy;
                    clear RS_index;
                    if sum(50<Pso.Gbest_value-Pso.Fit_var)~=0
                        dummy = 50 + Pso.Fit_var;
                        dummy(find(dummy>Pso.Gbest_value))=inf;
                        dummy(find(dummy<inf))=0;
                        [~,RS_index] = min (dummy + (Pso.vlength./max(Pso.vlength))+(Pso.heigth_var./max(Pso.heigth_var)));
                    else
                        [~,RS_index] = max(Pso.Gbest_value);
                    end
                    ROOT.solution_counter_T50  = ROOT.solution_counter_T50  + 1;
                    ROOT.solution_position_T50  = Pso.Gbest_position(RS_index,:);
                    ROOT.solution_value_T50  = Pso.Gbest_value(RS_index);
                    ROOT.survival_time_sequence_T50 (1,ceil((evals)/change_frequency))=0;
                else
                    ROOT.survival_time_sequence_T50 (1,ceil((evals)/change_frequency)) = ROOT.survival_time_sequence_T50 (1,((ceil((evals)/change_frequency))-1)) + 1;
                end
            end
            ROOT.solution_fitness_sequence_T40 (1,ceil((evals)/change_frequency)) = ROOT.solution_value_T40 ;
            %____________________________________________________
            ROOT.solution_fitness_sequence_T45 (1,ceil((evals)/change_frequency)) = ROOT.solution_value_T45 ;
            %____________________________________________________
            ROOT.solution_fitness_sequence_T50 (1,ceil((evals)/change_frequency)) = ROOT.solution_value_T50 ;
            %____________________________________________________
        end
        if (Pso.fe < Pso.max_fe)
            Pso.fe = evals;
%             Pso.FE.global_max(Pso.fe) = global_max-current_error;%
%             disp(['Run number = ' num2str(Pso.test) ', FE_number = ' num2str(Pso.fe)  ', Current_Error = ' num2str(current_error) ', Swarm# = ' num2str(Pso.TSN)]);
        else
            Pso.end = 1;
        end
    end
else
    Pso.end = 1;
end
end
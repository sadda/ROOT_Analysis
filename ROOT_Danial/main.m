close ;clear;clc;
format bank;
Pso.offline_err_FE=0;
for test=1 : 30
    %%%%%%%%%%%%%%%Initializing BEGIN
    clear Pso ROOT;
    Pso.test = test;
    global recent_change_PSO, global environment_change,global ROOT,global Last,global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType; %#ok<*TLEV,NUSED>
    init_parameters;
    init_peaks;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% ROOT_FE
    ROOT.solution_value_T40 = [];
    ROOT.solution_position_T40  = [];
    ROOT.solution_fitness_sequence_T40 = [];
    ROOT.survival_time_sequence_T40 = [];
    ROOT.solution_counter_T40 = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%
    ROOT.solution_value_T45 = [];
    ROOT.solution_position_T45  = [];
    ROOT.solution_fitness_sequence_T45 = [];
    ROOT.survival_time_sequence_T45 = [];
    ROOT.solution_counter_T45 = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ROOT.solution_value_T50 = [];
    ROOT.solution_position_T50  = [];
    ROOT.solution_fitness_sequence_T50 = [];
    ROOT.survival_time_sequence_T50 = [];
    ROOT.solution_counter_T50 = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pso.TSN = 0;% number of tracker swarms
    Pso.Q = 1;
    Pso.P = 1;
    Pso.ConvLimit = 1;
    Pso.TSps = 5;% tracker swarm population size
    Pso.excl_factor = 0.1;
    Pso.finder.K = 10;
    Pso.number = Pso.TSN * Pso.TSps;
    Pso.D=geno_size;
    Pso.max_fe=environment_change * change_frequency;
    Pso.FE.current_max = zeros(1 , Pso.max_fe);
    Pso.end = 0;
    Pso.Pbest_value = inf(Pso.number , 1);
    Pso.Pbest_position = zeros(Pso.number , Pso.D);
    Pso.fe = 0;
    Pso.old = [];
    Pso.Gbest_past_environment = [];
    Pso.velocity = zeros(Pso.number , Pso.D);
    Pso.change_accur = 0;
    Pso.excl = 5;
    Pso.init_vlength=1;
    Pso.vlength = [];
    Pso.itr = 0;
    Pso.fitness_function_bound = maxcoordinate;
    Pso.X = mincoordinate + ((maxcoordinate-mincoordinate) * (rand(Pso.number , Pso.D)));
    Pso.Pbest_position = Pso.X;
    Pso.particle_value = inf(Pso.number , 1);
    Pso.Gbest_value = inf(1,Pso.TSN);
    Pso.Gbest_position = inf(Pso.TSN , Pso.D);
    Pso.Gbest_past_value  = Pso.Gbest_value;
    Pso.x = 0.729843788;
    c1 = 2.05;
    c2 = 2.05;
    Pso.itr=0;
    Pso.heigth_var = [];
    Pso.Fit_var = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Finder
    Pso.finder.number = 10;
    Pso.finder.dimension = Pso.D;
    Pso.finder.X = mincoordinate + (rand(Pso.finder.number , Pso.finder.dimension) * (maxcoordinate-mincoordinate));
    Pso = fitness(Pso.finder.X , Pso);
    Pso.finder.value = Pso.result;
    Pso.finder.Pbest_position = Pso.finder.X;
    Pso.finder.Pbest_value = Pso.finder.value;
    [value, index] = max(Pso.finder.Pbest_value);
    Pso.finder.Gbest_value = value;
    Pso.finder.Gbest_position = Pso.finder.Pbest_position(index,:);
    Pso.finder.velocity = rands(Pso.finder.number , Pso.finder.dimension) * 10;
    Pso.finder.activate = 0;
    Pso.finder.Vmin = -1 * 0.3 * (maxcoordinate-mincoordinate);
    Pso.finder.Vmax = 0.3 * (maxcoordinate-mincoordinate);
    Pso.finder.Gbest_position_history = [];
    %%%%%%%%%%%%%%%Initializing END
%     figure(1); hold on;
%     [X,Y] = meshgrid(mincoordinate:1:maxcoordinate);
%     Z = (peak(1,end) - peak(1,end-1)*(sqrt((X-peak(1,1)).^2 + (Y-peak(1,2)).^2)));
    while(1==1)
%         drawnow;
%         hold off;
%         Z = (peak(1,end) - peak(1,end-1)*(sqrt((X-peak(1,1)).^2 + (Y-peak(1,2)).^2)));
%         kk = 2;
%         while kk <= number_of_peaks
%             ZI =(peak(kk,end) - peak(kk,end-1)*(sqrt((X-peak(kk,1)).^2 + (Y-peak(kk,2)).^2)));
%             Z = max(Z,ZI);
%             kk = kk + 1;
%         end
%         contour(X,Y,Z,20);
%         hold on;
%         for i=1 : Pso.TSN
%             plot(Pso.Gbest_position(i,1),Pso.Gbest_position(i,2),'*','markersize',10,'color','k');
%         end
%         plot(Pso.finder.Gbest_position(:,1),Pso.finder.Gbest_position(:,2),'p','markersize',10,'color','r','markerfacecolor','r');
%         plot(peak(maximum_peak,1),peak(maximum_peak,2),'o','markersize',5,'color','k','markerfacecolor','k');
        %% Test for change
        if (recent_change_PSO == 1)
            recent_change_PSO=0;
            disp(['Run number = ' num2str(Pso.test)]);
            %% re-evaluating robust solutions
            Pso = fitness(ROOT.solution_position_T40 , Pso);
            ROOT.solution_value_T40 = Pso.result;
            %____________________________________________
            Pso = fitness(ROOT.solution_position_T45 , Pso);
            ROOT.solution_value_T45 = Pso.result;
            %____________________________________________
            Pso = fitness(ROOT.solution_position_T50 , Pso);
            ROOT.solution_value_T50 = Pso.result;
            %____________________________________________    
            %% finder after change detection
            Pso = fitness(Pso.finder.Pbest_position , Pso);
            Pso.finder.Pbest_value = Pso.result;
            Pso = fitness(Pso.finder.X , Pso);
            Pso.finder.value = Pso.result;
            better = (Pso.finder.value > Pso.finder.Pbest_value);
            for j = 1: Pso.finder.number
                if(better(j))
                    Pso.finder.Pbest_value(j) = Pso.finder.value(j);
                    Pso.finder.Pbest_position(j , :) = Pso.finder.X(j , :);
                end
            end
            [Pso.finder.Gbest_value, index] = max(Pso.finder.Pbest_value);
            Pso.finder.Gbest_position = Pso.finder.Pbest_position(index , :);
            %% Shift Severiry
            Pso.last_vlength = [];
            for i=1 : Pso.TSN
                if Pso.old(i) >1
                    Pso.last_vlength(i) = Euc_dis(Pso.Gbest_past_environment(i,:) , Pso.Gbest_position(i,:));
                    Pso.vlength(i) = ((Pso.old(i)*Pso.vlength(i)) + Pso.last_vlength(i))/(Pso.old(i)+1);
                elseif Pso.old(i) == 1
                    Pso.vlength(i) = Euc_dis(Pso.Gbest_past_environment(i,:) , Pso.Gbest_position(i,:));
                    Pso.last_vlength(i) = Pso.vlength(i);
                elseif Pso.old(i) == 0
                    Pso.last_vlength(i) = 1;
                end
            end
            %% heigth_var
            for i=1 : Pso.TSN
                if Pso.old(i) > 1
                    Pso.heigth_var(i) = ((Pso.heigth_var(i)*Pso.old(i)) + abs(abs(Pso.Gbest_past_value(i)) - abs(Pso.Gbest_value(i))))/(Pso.old(i)+1);
                elseif Pso.old(i) == 1
                    Pso.heigth_var(i) = abs(abs(Pso.Gbest_past_value(i)) - abs(Pso.Gbest_value(i)));
                end
            end
            %% Fit_var
            Pso.Gbest_past_value = Pso.Gbest_value;
            Pso.Gbest_past_environment = Pso.Gbest_position;
            Pso = fitness(Pso.Gbest_position,Pso);
            dummy = Pso.result;
            for i=1 : Pso.TSN
                if Pso.old(i) > 1
                    Pso.Fit_var(i) =  ((Pso.Fit_var(i) * Pso.old(i)) + (abs(Pso.Gbest_past_value(i)) - abs(dummy(i)))) / (Pso.old(i)+1);
                elseif Pso.old(i) == 1
                    Pso.Fit_var(i) = abs(abs(Pso.Gbest_past_value(i)) - abs(dummy(i)));
                end
            end
            %%%%%%%%%%%%%%%%%%%
            for i=1 : Pso.TSN
                Pso.velocity((((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ) , :)= rands(Pso.TSps  , Pso.D) * (Pso.Q* Pso.vlength(i));
            end
            for i=1 : Pso.TSN
                Pso.X((((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ) , :)...
                    = repmat(Pso.Gbest_position(i,:), Pso.TSps , 1) + rands(Pso.TSps  , Pso.D) * (Pso.P * Pso.vlength(i));
                Pso.Pbest_position((((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ) , :)= Pso.X((((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ) , :);
                Pso = fitness(Pso.Pbest_position(( ((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ) , :) , Pso);
                Pso.Pbest_value((((i-1)*Pso.TSps + 1) : (i*Pso.TSps))) = Pso.result;
            end
            for i=1:Pso.TSN
                [val, id] = max(Pso.Pbest_value( ((i-1)*Pso.TSps + 1) : (i*Pso.TSps) ));
                id = id + ((i-1)*Pso.TSps);
                Pso.Gbest_value(i) = val;
                Pso.Gbest_position(i,:) = Pso.Pbest_position(id,:);
            end
            Pso.change_accur = 0;
            Pso.old = Pso.old+ones(size(Pso.old));
        end
        Pso.itr = Pso.itr + 1;
        %% Finder Movement begin
        Gbest_temp = repmat(Pso.finder.Gbest_position , Pso.finder.number , 1);
        Pso.finder.velocity = Pso.x * ((Pso.finder.velocity) + (c1*rand(Pso.finder.number , Pso.finder.dimension).*(Pso.finder.Pbest_position - Pso.finder.X))+ (c2*rand(Pso.finder.number , Pso.finder.dimension).*(Gbest_temp - Pso.finder.X)));
        Pso.finder.velocity(find( Pso.finder.velocity > Pso.finder.Vmax )) = Pso.finder.Vmax; %#ok<*FNDSB>
        Pso.finder.velocity(find( Pso.finder.velocity < Pso.finder.Vmin )) = Pso.finder.Vmin;
        Pso.finder.X = Pso.finder.X + Pso.finder.velocity;
        Pso = fitness(Pso.finder.X , Pso);
        Pso.finder.value = Pso.result;
        better = (Pso.finder.value > Pso.finder.Pbest_value);
        for j = 1: Pso.finder.number
            if(better(j))
                Pso.finder.Pbest_value(j) = Pso.finder.value(j);
                Pso.finder.Pbest_position(j , :) = Pso.finder.X(j , :);
            end
        end
        [Pso.finder.Gbest_value, index] = max(Pso.finder.Pbest_value);
        Pso.finder.Gbest_position = Pso.finder.Pbest_position(index , :);
        Pso.finder.Gbest_val(Pso.itr) = Pso.finder.Gbest_value;
        Pso.finder.Gbest_position_history = [Pso.finder.Gbest_position_history ; Pso.finder.Gbest_position;];
        [row,~] = size(Pso.finder.Gbest_position_history);
        if row > Pso.finder.K
            Pso.finder.Gbest_position_history = removerows(Pso.finder.Gbest_position_history, 1);
        end
        %% Finder Exclusion BEGIN
        for i=1:Pso.TSN
            if (Euc_dis(Pso.Gbest_position(i,:) , Pso.finder.Gbest_position) <= Pso.excl)
                if Pso.finder.Gbest_value > Pso.Gbest_value(i)
                    Pso.Gbest_value(i) = Pso.finder.Gbest_value;
                    Pso.Gbest_position(i,:) = Pso.finder.Gbest_position;
                end
                Pso.finder.X = mincoordinate + (rand(Pso.finder.number , Pso.finder.dimension) * (maxcoordinate-mincoordinate));
                Pso = fitness(Pso.finder.X , Pso);
                Pso.finder.value = Pso.result;
                Pso.finder.Pbest_position = Pso.finder.X;
                Pso.finder.Pbest_value = Pso.finder.value;
                [value, index] = min(Pso.finder.Pbest_value);
                Pso.finder.Gbest_value = value;
                Pso.finder.Gbest_position = Pso.finder.Pbest_position(index,:);
                Pso.finder.velocity = rands(Pso.finder.number , Pso.finder.dimension);
                break;
            end
        end
        %%  finder Convergence BEGIN
        if Pso.itr>Pso.finder.K
            if Euc_dis(Pso.finder.Gbest_position_history(1 , :) , Pso.finder.Gbest_position_history(Pso.finder.K , :)) < Pso.ConvLimit 
                Pso.finder.activate = 1;
            end
        end
        %% Execute activation BEGIN
        if Pso.finder.activate == 1
            [~, index_list] = sort(Pso.finder.Pbest_value,'descend');
            %%%creat new swarm
            Pso.TSN = Pso.TSN + 1;
            Pso.vlength = [Pso.vlength ; Pso.init_vlength];
            Pso.heigth_var = [Pso.heigth_var ; inf];
            Pso.Fit_var = [Pso.Fit_var ; inf];
            Pso.old = [Pso.old ; 0];
            Pso.Gbest_value = [Pso.Gbest_value ; Pso.finder.Gbest_value];
            Pso.Gbest_position = [Pso.Gbest_position ; Pso.finder.Gbest_position];
            Pso.Gbest_past_environment = [Pso.Gbest_past_environment ; Pso.finder.Gbest_position];
            Pso.X = [Pso.X ; Pso.finder.X(index_list(1:Pso.TSps) , :);];
            Pso.particle_value = [Pso.particle_value ; Pso.finder.value(index_list(1:Pso.TSps))];
            Pso.velocity = [Pso.velocity ; Pso.finder.velocity(index_list(1:Pso.TSps) , :)];
            Pso.Pbest_value = [Pso.Pbest_value ; Pso.finder.Pbest_value(index_list(1:Pso.TSps))];
            Pso.Pbest_position = [Pso.Pbest_position ; Pso.finder.Pbest_position(index_list(1:Pso.TSps),:)];
            %%%Reset finder
            Pso.finder.X = mincoordinate + (rand(Pso.finder.number , Pso.finder.dimension) * (maxcoordinate-mincoordinate));
            Pso = fitness(Pso.finder.X , Pso);
            Pso.finder.value = Pso.result;
            Pso.finder.Pbest_position = Pso.finder.X;
            Pso.finder.Pbest_value = Pso.finder.value;
            [value, index] = max(Pso.finder.Pbest_value);
            Pso.finder.Gbest_value = value;
            Pso.finder.Gbest_position = Pso.finder.Pbest_position(index,:);
            Pso.finder.velocity = rands(Pso.finder.number , Pso.finder.dimension) * 20;
            Pso.finder.activate = 0;
            Pso.excl = Pso.excl_factor * ((maxcoordinate-mincoordinate) / ((Pso.TSN) ^ (1 / Pso.D)));   % r_excl = 0.5*d_boa
        end
        %% Update  particles BEGIN
        for i=1 : Pso.TSN
            if recent_change_PSO==1
                break;
            end
            for j=1:Pso.TSps
                Pso.velocity((i-1)*Pso.TSps + j , :) = ...
                    Pso.x * ((Pso.velocity((i-1)*Pso.TSps + j , :) + ...
                    (c1*rand(1 , Pso.D).*(Pso.Pbest_position((i-1)*Pso.TSps + j , :) - Pso.X((i-1)*Pso.TSps + j , :))) + (c2*rand(1 , Pso.D).*(Pso.Gbest_position(i,:) - Pso.X((i-1)*Pso.TSps + j , :)))));
                if Pso.old(i)>1
                    for k=1:Pso.D
                        if Pso.velocity((i-1)*Pso.TSps + j , k) > Pso.vlength(i)
                            Pso.velocity((i-1)*Pso.TSps + j , k) = Pso.vlength(i);
                        elseif Pso.velocity((i-1)*Pso.TSps + j , k) < (-1 * Pso.vlength(i))
                            Pso.velocity((i-1)*Pso.TSps + j , k) = -1 * Pso.vlength(i);
                        end
                    end
                end
                Pso.X((i-1)*Pso.TSps + j,:) = Pso.X((i-1)*Pso.TSps + j , :) + Pso.velocity((i-1)*Pso.TSps + j , :);
                Pso = fitness(Pso.X((((i-1) * Pso.TSps) + j)  , :) , Pso);
                Pso.particle_value((((i-1) * Pso.TSps) + j) ) = Pso.result;
            end
            for j=1 : Pso.TSps
                if (Pso.particle_value(((i-1) * Pso.TSps) + j) > Pso.Pbest_value(((i-1) * Pso.TSps) + j))
                    Pso.Pbest_value(((i-1) * Pso.TSps) + j) = Pso.particle_value(((i-1) * Pso.TSps) + j);
                    Pso.Pbest_position((((i-1) * Pso.TSps) + j) , :) = Pso.X((((i-1) * Pso.TSps) + j) , :);
                end
            end
            [val, id] = max(Pso.Pbest_value((((i-1) * Pso.TSps) + 1) : (i * Pso.TSps)));
            Pso.Gbest_value(i) = val;
            Pso.Gbest_position(i,:) = Pso.Pbest_position((((i-1) * Pso.TSps) + id), :);
        end
        %% trackers exclusion BEGIN
        flag = 0;
        for i=1:Pso.TSN
            for j=1:Pso.TSN
                if flag==0
                    if (i ~= j) && (Euc_dis(Pso.Gbest_position(i,:) , Pso.Gbest_position(j,:)) <= Pso.excl)
                        if Pso.Gbest_value(i)> Pso.Gbest_value(j)
                            kill = j;
                        else
                            kill = i;
                        end
                        if Pso.old(j) < Pso.old(i) && kill==i
                            Pso.vlength(j) = Pso.vlength(i);
                            Pso.heigth_var(j) = Pso.heigth_var(i);
                            Pso.Fit_var(j) = Pso.heigth_var(i);
                        elseif Pso.old(j) > Pso.old(i) && kill==j
                            Pso.vlength(i) = Pso.vlength(j);
                            Pso.heigth_var(i) = Pso.heigth_var(j);
                            Pso.Fit_var(i) = Pso.heigth_var(j);
                        end
                        Pso.TSN = Pso.TSN - 1;
                        Pso.vlength = removerows(Pso.vlength , kill);
                        Pso.old = removerows(Pso.old , kill);
                        Pso.Gbest_value = removerows(Pso.Gbest_value , kill);
                        Pso.Gbest_position = removerows(Pso.Gbest_position , kill);
                        Pso.Gbest_past_environment = removerows(Pso.Gbest_past_environment , kill);
                        Pso.X = removerows(Pso.X , (((kill-1)*Pso.TSps)+1 : (kill*Pso.TSps)));
                        Pso.particle_value = removerows(Pso.particle_value , (((kill-1)*Pso.TSps)+1 : (kill*Pso.TSps)));
                        Pso.velocity = removerows(Pso.velocity , (((kill-1)*Pso.TSps)+1 : (kill*Pso.TSps)));
                        Pso.Pbest_value = removerows(Pso.Pbest_value , (((kill-1)*Pso.TSps)+1 : (kill*Pso.TSps)));
                        Pso.Pbest_position = removerows(Pso.Pbest_position , (((kill-1)*Pso.TSps)+1 : (kill*Pso.TSps)));
                        Pso.heigth_var = removerows(Pso.heigth_var , kill);
                        Pso.Fit_var = removerows(Pso.Fit_var , kill);
                        flag = 1;
                        Pso.excl = Pso.excl_factor * ((maxcoordinate-mincoordinate) / ((Pso.TSN) ^ (1 / Pso.D)));   % r_excl = 0.5*d_boa
                    end
                end
            end
        end
        %%
        if (Pso.end == 1)
            XXX.ROOT_value_in_each_environment_T40(test,:) = ROOT.solution_fitness_sequence_T40;
            XXX.ROOT_solution_num_T40(test) = ROOT.solution_counter_T40;
            XXX.ROOT_survival_time_sequence_T40(test,:) = ROOT.survival_time_sequence_T40;
            %________________________________________________
            XXX.ROOT_value_in_each_environment_T45(test,:) = ROOT.solution_fitness_sequence_T45;
            XXX.ROOT_solution_num_T45(test) = ROOT.solution_counter_T45;
            XXX.ROOT_survival_time_sequence_T45(test,:) = ROOT.survival_time_sequence_T45;
            %________________________________________________
            XXX.ROOT_value_in_each_environment_T50(test,:) = ROOT.solution_fitness_sequence_T50;
            XXX.ROOT_solution_num_T50(test) = ROOT.solution_counter_T50;
            XXX.ROOT_survival_time_sequence_T50(test,:) = ROOT.survival_time_sequence_T50;
            %________________________________________________
            break;
        end
    end
    %%
end
%%
Offline_error = mean(mean(Last.error));
%________________________________________________________
ROOT_fitness_value_T40 = mean(mean(XXX.ROOT_value_in_each_environment_T40'));
ROOT_fitness_value_T40_Std_Err = std(mean(XXX.ROOT_value_in_each_environment_T40'))/sqrt(test);
ROOT_solution_number_T40 = mean(XXX.ROOT_solution_num_T40);
ROOT_solution_life_cycle_T40 = environment_change/ROOT_solution_number_T40;
ROOT_solution_life_cycle_T40_Std_Err = std(environment_change./XXX.ROOT_solution_num_T40)/sqrt(test);
ROOT_survival_time_T40 = mean(mean(XXX.ROOT_survival_time_sequence_T40'));
ROOT_survival_time_T40_Std_Err = std(mean(XXX.ROOT_survival_time_sequence_T40'))/sqrt(test);
%________________________________________________________
ROOT_fitness_value_T45 = mean(mean(XXX.ROOT_value_in_each_environment_T45'));
ROOT_solution_number_T45 = mean(XXX.ROOT_solution_num_T45);
ROOT_solution_life_cycle_T45 = environment_change/ROOT_solution_number_T45;
ROOT_solution_life_cycle_T45_Std_Err = std(environment_change./XXX.ROOT_solution_num_T45)/sqrt(test);
ROOT_survival_time_T45 = mean(mean(XXX.ROOT_survival_time_sequence_T45'));
ROOT_survival_time_T45_Std_Err = std(mean(XXX.ROOT_survival_time_sequence_T45'))/sqrt(test);
%________________________________________________________
ROOT_fitness_value_T50 = mean(mean(XXX.ROOT_value_in_each_environment_T50'));
ROOT_fitness_value_T50_Std_Err = std(mean(XXX.ROOT_value_in_each_environment_T50'))/sqrt(test);
ROOT_solution_number_T50 = mean(XXX.ROOT_solution_num_T50);
ROOT_solution_life_cycle_T50 = environment_change/ROOT_solution_number_T50;
ROOT_solution_life_cycle_T50_Std_Err = std(environment_change./XXX.ROOT_solution_num_T50)/sqrt(test);
ROOT_survival_time_T50 = mean(mean(XXX.ROOT_survival_time_sequence_T50'));
ROOT_survival_time_T50_Std_Err = std(mean(XXX.ROOT_survival_time_sequence_T50'))/sqrt(test);
%________________________________________________________
[ROOT_survival_time_T40,ROOT_survival_time_T40_Std_Err;ROOT_survival_time_T45,ROOT_survival_time_T45_Std_Err;ROOT_survival_time_T50,ROOT_survival_time_T50_Std_Err]

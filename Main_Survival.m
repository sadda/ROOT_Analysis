clear all;
close all;

addpath(genpath('../ROOT-Benchmark'));

opts = Initialize_Options('Default1');

S = 10000;
m = 5;

c = opts.x_min + (opts.x_max - opts.x_min)*rand(S,m);
h = opts.h_min + (opts.h_max - opts.h_min)*rand(S,m);
w = opts.w_min + (opts.w_max - opts.w_min)*rand(S,m);

% delta = opts.h_min;
delta = 50;

[~, i_tmo] = max(h, [], 2);
[~, i_rob] = max((h-delta)./w, [], 2);

n_max = 100;

x = opts.x_min + (opts.x_max - opts.x_min)*rand(S,n_max);


f_aux = zeros(S,n_max,m);
for n=1:n_max
    f_aux(:,n,:) = h - w.*abs(c - x(:,n));
end

[f_max, i_max] = max(f_aux, [], 3);




i_our = zeros(S,n_max);
for n=1:n_max
    
    [~, i_row] = max(f_max(:,1:n), [], 2);
    i_our(:,n) = Select_Rows(i_max, i_row);
    
end




Surv_Time = @(ii, h, w, s, delta) max(ceil((Select_Rows(h, ii) - delta) ./ (s*Select_Rows(w, ii))), 0).^2;

s = 3;


time_tmo = Surv_Time(i_tmo, h, w, s, delta);
time_rob = Surv_Time(i_rob, h, w, s, delta);
time_our = zeros(S,n_max);
for n=1:n_max
    time_our(:,n) = Surv_Time(i_our(:,n), h, w, s, delta);
end



fig = figure();
plot(1:n_max, repmat(mean(i_rob == i_tmo), n_max, 1));
hold on;
plot(1:n_max, mean(i_our == i_rob));
plot(1:n_max, mean(i_our == i_tmo));
legend({'Robust = TMO', 'Ours = Robust', 'Ours = TMO'}, 'location', 'northwest');
xlabel('Number of function evaluations');
ylabel('Percentage of equality');
title(sprintf('Threshold = %d', delta));
saveas(fig, sprintf('Results_Delta=%d_1.jpg', delta));


fig = figure();
plot(1:n_max, repmat(mean(time_rob), n_max, 1));
hold on;
plot(1:n_max, mean(time_our));
plot(1:n_max, repmat(mean(time_tmo), n_max, 1));
legend({'Robust', 'Ours', 'TMO'});
xlabel('Number of function evaluations');
ylabel('Survival time');
title(sprintf('Threshold = %d', delta));
saveas(fig, sprintf('Results_Delta=%d_2.jpg', delta));





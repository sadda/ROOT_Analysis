clear all;
close all;

addpath('../');


seeds = [6 8 8];

for i_seed = 1:length(seeds)
    
    rng(seeds(i_seed));
    
    m  = 5;
    n1 = 11;
    n2 = 101;
    
    x_min = 0;
    x_max = 10;
    
    
    c = x_min + (x_max-x_min)*rand(m,2);
    h = 30+(70-30)*rand(m,1);
    w = 1+(13-1)*rand(m,1);
    f = @(x) Compute_F(x, h, w, c, 'bench1');
    
    if i_seed <= 2
        [X1, Y1] = meshgrid(linspace(x_min, x_max, n1), linspace(x_min, x_max, n1));
        X1 = X1(:);
        Y1 = Y1(:);
    else
        rng(seeds(i_seed)+4);
        
        X1 = x_min + (x_max-x_min)*rand(n1^2,1);
        Y1 = x_min + (x_max-x_min)*rand(n1^2,1);
    end
    
    Z1 = f([X1,Y1]);
    
    [X2, Y2] = meshgrid(linspace(x_min, x_max, n2), linspace(x_min, x_max, n2));
    X2 = X2(:);
    Y2 = Y2(:);
    Z2 = f([X2,Y2]);
    
    Z_max = max(Z2);
    
    [~,i_max1] = max(Z1);
    [~,i_max2] = max(Z2);
    
    dist     = vecnorm([X1(i_max1), Y1(i_max1)] - c, 2, 2);
    [~,opt1] = min(dist);
    [~,opt2] = max(h);
    
    
    
    fig = figure();
    hold on;
    
    gd = delaunay(X2, Y2);
    
    trisurf(gd, X2, Y2, Z2);
    view(2);
    shading interp;
    colormap jet;
    colorbar;
    axis off;
    scatter3(c(opt2,1), c(opt2,2), Z_max, 80, 'k', '*');
    
    file_name = sprintf('Res_%d_%d.jpg', i_seed, 1);
    saveas(fig, file_name);
    CropImage(file_name, file_name);
    
    scatter3(X1, Y1, repmat(Z_max, size(X1)), 10, 'k', 'filled');
    
    file_name = sprintf('Res_%d_%d.jpg', i_seed, 2);
    saveas(fig, file_name);
    CropImage(file_name, file_name);
    
    scatter3(X1(i_max1), Y1(i_max1), Z_max, 40, 'k', 'filled');
    
    file_name = sprintf('Res_%d_%d.jpg', i_seed, 3);
    saveas(fig, file_name);
    CropImage(file_name, file_name);
    
    scatter3(c(opt1,1), c(opt1,2), Z_max, 40, 'k', 'filled');
    
    file_name = sprintf('Res_%d_%d.jpg', i_seed, 4);
    saveas(fig, file_name);
    CropImage(file_name, file_name);
    
end


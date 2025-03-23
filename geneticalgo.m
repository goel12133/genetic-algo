function [I_opt, best_I, best_cost] = tokamak_optimization()
clear; clc; close all;

mu0 = 4*pi*1e-7;
R0 = 2.5;
a = 1.0;
n_coils = 12;
n_poloidal = 6;

critical_points = [
    R0 + 0.9*a, 0.5*a;
    R0 + 0.9*a, -0.5*a;
    R0 - 0.9*a, 0.5*a;
    R0 - 0.9*a, -0.5*a;
    R0, 0.0;
];

[tor_coil_pos, tor_coil_dir, pol_coil_pos, pol_coil_dir] = ...
    generate_coils(R0, a, n_coils, n_poloidal);

n_poloidal = size(pol_coil_pos, 1);

I_toroidal = 1e6 * ones(n_coils, 1);
I_poloidal = 5e5 * ones(n_poloidal, 1);

[R_grid,Z_grid] = meshgrid(linspace(R0-a,R0+a,50),linspace(-a,a,50));

target.B_tor_target = 5.0 * ones(size(R_grid));
target.B_pol_target = 0.5 * ones(size(R_grid));

I0 = [I_toroidal; I_poloidal];
lb = 0.2*I0;
ub = 5.0*I0;

options = optimoptions('ga',...
    'Display', 'iter',...
    'PopulationSize', 100,...
    'MaxGenerations', 200,...
    'FunctionTolerance', 1e-6,...
    'CrossoverFraction', 0.8,...
    'MutationFcn', {@mutationadaptfeasible, 0.2},...
    'UseVectorized', false,...
    'OutputFcn', @outfun);

[I_opt, fval] = ga(@(I)optimization_cost(I, target, critical_points, R_grid, Z_grid,...
                        tor_coil_pos, tor_coil_dir,...
                        pol_coil_pos, pol_coil_dir,...
                        n_coils, n_poloidal, mu0),...
                        length(I0), [], [], [], [], lb, ub, [], options);

best_I = I_opt;
best_cost = fval;

plot_results(I_opt, R_grid, Z_grid, critical_points,...
            tor_coil_pos, tor_coil_dir,...
            pol_coil_pos, pol_coil_dir,...
            n_coils, n_poloidal, mu0);

end

function [tor_pos, tor_dir, pol_pos, pol_dir] = generate_coils(R0, a, n_coils, n_poloidal)
    theta = linspace(0, 2*pi, n_coils+1);
    theta(end) = [];
    tor_pos = [R0*cos(theta'), R0*sin(theta'), zeros(n_coils,1)];
    tor_dir = repmat([0 0 1], n_coils, 1);
    
    z_pos = linspace(-a, a, n_poloidal);
    
    pol_pos_outer = [(R0 + a)*ones(n_poloidal,1), zeros(n_poloidal,1), z_pos'];
    pol_pos_inner = [(R0 - a)*ones(n_poloidal,1), zeros(n_poloidal,1), z_pos'];
    pol_pos_mid = [R0*ones(n_poloidal,1), zeros(n_poloidal,1), z_pos'];
    
    pol_pos = [pol_pos_outer; pol_pos_inner; pol_pos_mid];
    pol_dir = repmat([0 1 0], size(pol_pos,1), 1);
end

function B = compute_B_at_points(points, coil_pos, coil_dir, currents, mu0)
    B = zeros(size(points));
    for c = 1:size(coil_pos, 1)
        dl = coil_dir(c,:);
        dl_expanded = repmat(dl, size(points,1), 1);
        r_vec = points - coil_pos(c,:);
        r_mag = vecnorm(r_vec, 2, 2);
        cross_prod = cross(dl_expanded, r_vec, 2);
        dB = (mu0*currents(c)/(4*pi)) * cross_prod ./ r_mag.^3;
        B = B + dB;
    end
end

function B_fields = calculate_fields(I, points, coil_pos, coil_dir, n_coils, mu0)
    if length(I) ~= size(coil_pos, 1)
        error('Currents/coils mismatch: %d currents vs %d coils', length(I), size(coil_pos,1));
    end
    
    B = compute_B_at_points(points, coil_pos, coil_dir, I, mu0);
    B_fields = vecnorm(B, 2, 2);
end

function cost = optimization_cost(I, target, critical_points, R_grid, Z_grid,...
                                tor_pos, tor_dir, pol_pos, pol_dir,...
                                n_coils, n_poloidal, mu0)
    I_tor = I(1:n_coils);
    I_pol = I(n_coils+1:end);
    
    grid_points = [R_grid(:), Z_grid(:), zeros(numel(R_grid),1)];
    B_tor = calculate_fields(I_tor, grid_points, tor_pos, tor_dir, n_coils, mu0);
    B_pol = calculate_fields(I_pol, grid_points, pol_pos, pol_dir, n_poloidal, mu0);
    
    crit_points_3D = [critical_points, zeros(size(critical_points,1),1)];
    B_pol_crit = calculate_fields(I_pol, crit_points_3D, pol_pos, pol_dir, n_poloidal, mu0);
    
    cost_tor = mean((B_tor - target.B_tor_target(:)).^2);
    cost_pol = mean((B_pol - target.B_pol_target(:)).^2);
    
    cost_walls = sum(max(0, 1.0 - B_pol_crit(1:4)).^2);
    
    cost_center = min((B_pol_crit(5) - 0.2)^2, 0.1);
    
    cost = cost_tor + cost_pol + 10000 * cost_walls + 500 * cost_center;
end

function plot_results(I_opt, R_grid, Z_grid, critical_points,...
                    tor_pos, tor_dir, pol_pos, pol_dir,...
                    n_coils, n_poloidal, mu0)
    grid_points = [R_grid(:), Z_grid(:), zeros(numel(R_grid),1)];
    B_tor = calculate_fields(I_opt(1:n_coils), grid_points, tor_pos, tor_dir, n_coils, mu0);
    B_pol = calculate_fields(I_opt(n_coils+1:end), grid_points, pol_pos, pol_dir, n_poloidal, mu0);
    
    B_tor = reshape(B_tor, size(R_grid));
    B_pol = reshape(B_pol, size(R_grid));
    
    crit_points_3D = [critical_points, zeros(size(critical_points,1),1)];
    B_pol_crit = calculate_fields(I_opt(n_coils+1:end), crit_points_3D, pol_pos, pol_dir, n_poloidal, mu0);
    
    figure('Position', [100 100 1200 800]);
    
    subplot(2,1,1);
    contourf(R_grid, Z_grid, B_tor);
    title('Optimized Toroidal Field');
    colorbar;
    axis equal;
    
    subplot(2,1,2);
    contourf(R_grid, Z_grid, B_pol);
    hold on;
    plot(critical_points(1:4,1), critical_points(1:4,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(critical_points(5,1), critical_points(5,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    title('Optimized Poloidal Field with Critical Points');
    legend('Wall points', 'Plasma center', 'Location', 'northeast');
    colorbar;
    axis equal;
    
    fprintf('\nPoloidal field at critical points:\n');
    fprintf('Outer wall upper: %.2f T\n', B_pol_crit(1));
    fprintf('Outer wall lower: %.2f T\n', B_pol_crit(2));
    fprintf('Inner wall upper: %.2f T\n', B_pol_crit(3));
    fprintf('Inner wall lower: %.2f T\n', B_pol_crit(4));
    fprintf('Plasma center:    %.2f T\n\n', B_pol_crit(5));
end

function [stop, options, optchanged] = outfun(optimValues, state, options)
    persistent best_I;  
    persistent best_cost;  
    persistent n_coils;  

    stop = false;       
    optchanged = false; 

    if isempty(n_coils)
        n_coils = evalin('base', 'n_coils');
    end

    if strcmp(state, 'iter')
        if isempty(best_I) || optimValues.fval < best_cost
            best_I = optimValues.bestx;  
            best_cost = optimValues.fval;  
        end
        
        fprintf('Generation %d: Best Cost = %.4f\n', optimValues.generation, best_cost);
    end

    if strcmp(state, 'done')
        fprintf('\nOptimization complete. Best I value:\n');
        fprintf('Toroidal currents: [');
        fprintf('%.2e ', best_I(1:n_coils));
        fprintf(']\nPoloidal currents: [');
        fprintf('%.2e ', best_I(n_coils+1:end));
        fprintf(']\nBest Cost: %.4f\n', best_cost);
    end
end

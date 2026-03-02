%% Task 2: Calculate gamma1 and gamma2 for Minimum Phase (Numerical Sweep)
resolution = 0.02;
gamma1_vec = 0.01:resolution:0.99; 
gamma2_vec = 0.01:resolution:0.99;

[G1_grid, G2_grid] = meshgrid(gamma1_vec, gamma2_vec);
MinPhaseRegion = zeros(size(G1_grid));

% Loop through all combinations
for i = 1:size(G1_grid, 1)
    for j = 1:size(G1_grid, 2)
        g1_val = G1_grid(i, j);
        g2_val = G2_grid(i, j);
        
        % Define G(s) for the current gamma pair
        G11_temp = (g1_val * c1_P) / (T1_P*s + 1);
        G12_temp = ((1 - g2_val) * c1_P) / ((T1_P*s + 1) * (T3_P*s + 1));
        
        % Using the denominator from the PDF: (T2*s + 1)*(T1*s + 1)
        G21_temp = ((1 - g1_val) * c2_P) / ((T2_P*s + 1) * (T1_P*s + 1)); 
        G22_temp = (g2_val * c2_P) / (T2_P*s + 1);
        
        G_current = [G11_temp, G12_temp; 
                     G21_temp, G22_temp];
        
        % Calculate transmission zeros
        z = tzero(G_current);
        
        % Check if all zeros are in the Left Half-Plane (LHP)
        if isempty(z)
            MinPhaseRegion(i, j) = 1; 
        elseif max(real(z)) < 0
            % Minimum phase (all zeros have negative real parts)
            MinPhaseRegion(i, j) = 1;
        else
            % Non-minimum phase (at least one zero has positive real part)
            MinPhaseRegion(i, j) = 0;
        end
    end
end

%% Plotting the Results
figure;
hold on;

% Plot the regions
% Green for Minimum Phase (1), Red for Non-Minimum Phase (0)
pcolor(G1_grid, G2_grid, MinPhaseRegion);
shading flat;
colormap([1 0.8 0.8; 0.8 1 0.8]); 

% Draw the theoretical boundary line: gamma1 + gamma2 = 1
plot([0 1], [1 0], 'k--', 'LineWidth', 2);

% Plot the P_minus operating point from the lab manual 
plot(0.70, 0.60, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8); % 
text(0.72, 0.62, 'P_- (Min Phase)', 'FontSize', 11, 'FontWeight', 'bold');

% Plot the P_plus operating point from the lab manual
plot(0.43, 0.34, 'rx', 'MarkerSize', 10, 'LineWidth', 3); % 
text(0.45, 0.36, 'P_+ (Non-Min Phase)', 'FontSize', 11, 'FontWeight', 'bold');

% Formatting
xlabel('\gamma_1');
ylabel('\gamma_2');
title('System Phase Behavior based on Flow Split Ratios');
xlim([0 1]);
ylim([0 1]);
grid on;
box on;

%% System Parameters (Nominal)
% [cm^2] Cross-section of tank 1-4 
A1 = 28;    
A3 = 28;    
A2 = 32;    
A4 = 32;   

% [cm^2] Cross-section of the outlet hole 1-4 
a1 = 0.071; 
a3 = 0.071; 
a2 = 0.057; 
a4 = 0.057; 

kc = 0.5;   % [V/cm] Sensor static gain 
g = 982;    % [cm/s^2] Acceleration of gravity 

% Tank water level initial values, tank 1-4
h1_init = 0;
h2_init = 0;
h3_init = 0;
h4_init = 0;

%% Operating Point Selection
% Select P_M variables for the linear model generation
h10_M = 12.3;
h20_M = 12.8;
h30_M = 1.6;
h40_M = 1.4;

h10_P = 12.3;
h20_P = 12.8;
h30_P = 1.6;
h40_P = 1.4;

k1_M = 3.33;
k2_M = 3.35;
gamma1_M = 0.7;
gamma2_M = 0.6;

k1_P = 3.14;
k2_P = 3.29;
gamma1_P = 0.43;
gamma2_P = 0.34;
%% Linear Model Generation (Task 1)
% Calculate time constants T_i 
T1_M = (A1/a1) * sqrt(2 * h10_M / g); 
T2_M = (A2/a2) * sqrt(2 * h20_M / g);
T3_M = (A3/a3) * sqrt(2 * h30_M / g);
T4_M = (A4/a4) * sqrt(2 * h40_M / g); 

T1_P = (A1/a1) * sqrt(2 * h10_P / g); 
T2_P = (A2/a2) * sqrt(2 * h20_P / g);
T3_P = (A3/a3) * sqrt(2 * h30_P / g);
T4_P = (A4/a4) * sqrt(2 * h40_P / g); 

c1_M = (T1_M * k1_M * kc) / A1;
c2_M = (T2_M * k2_M * kc) / A2;

c1_P = (T1_P * k1_P * kc) / A1;
c2_P = (T2_P * k2_P * kc) / A2;

s = tf('s');

% Define Transfer Function Matrix G(s) based on Eq. 7 in the lab manual 
G11_M = (gamma1_M * c1_M) / (T1_M*s + 1);
G12_M = ((1 - gamma2_M) * c1_M) / ((T1_M*s + 1) * (T3_M*s + 1)); 
G21_M = ((1 - gamma1_M) * c2_M) / ((T2_M*s + 1) * (T4_M*s + 1)); 
G22_M = (gamma2_M * c2_M) / (T2_M*s + 1);

G11_P = (gamma1_P * c1_P) / (T1_P*s + 1);
G12_P = ((1 - gamma2_P) * c1_P) / ((T1_P*s + 1) * (T3_P*s + 1)); 
G21_P = ((1 - gamma1_P) * c2_P) / ((T2_P*s + 1) * (T4_P*s + 1)); 
G22_P = (gamma2_P * c2_P) / (T2_P*s + 1);

% Combine into full MIMO transfer function
G_M = [G11_M, G12_M;
       G21_M, G22_M];

G_P = [G11_P , G12_P; 
       G21_P, G22_P];

% Calculate Poles and Zeros
sys_poles_M = pole(G_M);
sys_zeros_M = tzero(G_M)

sys_poles_P = pole(G_P);
sys_zeros_P = tzero(G_P)



% Task 2:
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


%% Conclusion -> if sum(gamma1, gamma2) > 1 -> G is a minimum phase system.
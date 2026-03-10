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
Ts = 0.01;
Delay = 0;

% Operating Point Selection
h10_M = 12.3;
h20_M = 12.8;
h30_M = 1.6;
h40_M = 1.4;
u_10_M = 3;
u_20_M = 3;

k1_M = 3.33;
k2_M = 3.35;
gamma1_M = 0.7;
gamma2_M = 0.6;

h10_P = 12.6;
h20_P = 13.0;
h30_P = 4.8;
h40_P = 4.9;
u_10_P = 3.15;
u_20_P = 3.15;

k1_P = 3.14;
k2_P = 3.29;
gamma1_P = 0.43;
gamma2_P = 0.34;

gamma1 = gamma1_M;
gamma2 = gamma2_M;

k1 = k1_M;
k2 = k2_M;

%% Tank water level initial values, tank 1-4
h1_init = h10_M;
h2_init = h20_M;
h3_init = h30_M;
h4_init = h40_M;

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
       G21_M, G22_M]

G_P = [G11_P , G12_P; 
       G21_P, G22_P];

% Calculate Poles and Zeros
sys_poles_M = pole(G_M);
sys_zeros_M = tzero(G_M);

sys_poles_P = pole(G_P);
sys_zeros_P = tzero(G_P);


%% Task 2:
%Task_2_plots
% Conclusion -> if sum(gamma1, gamma2) > 1 -> G is a minimum phase system.

%% Task 3: IMC
G11_num = cell2mat(G11_M.Numerator);
G11_den = cell2mat(G11_M.Denominator);
G11_inv = inv(G11_M);
lambda = 3;
Q11 = G11_inv/(lambda*s + 1)^2;
Q11_num = cell2mat(Q11.Numerator);
Q11_den = cell2mat(Q11.Denominator);

G22_num = cell2mat(G22_M.Numerator);
G22_den = cell2mat(G22_M.Denominator);
G22_inv = inv(G22_M);
lambda = 3;
Q22 = G22_inv/(lambda*s + 1)^2;
Q22_num = cell2mat(Q22.Numerator);
Q22_den = cell2mat(Q22.Denominator);

F_r_tilde_num = 1;
F_r_tilde_den = 1;

%% Task 4: LQG

% Calculate State Jacobian Matrix (A)
A = zeros(4, 4);

sq1 = sqrt(g / (2 * h10_M));
sq2 = sqrt(g / (2 * h20_M));
sq3 = sqrt(g / (2 * h30_M));
sq4 = sqrt(g / (2 * h40_M));

A(1,1) = -(a1 / A1) * sq1;
A(1,3) =  (a3 / A3) * sq3;
A(2,2) = -(a2 / A2) * sq2;
A(2,4) =  (a4 / A4) * sq4;
A(3,3) = -(a3 / A3) * sq3;
A(4,4) = -(a4 / A4) * sq4;

% Calculate Input Jacobian Matrix (B)
B = zeros(4, 2);

B(1,1) = (gamma1 * k1) / A1;
B(2,2) = (gamma2 * k2) / A2;
B(3,2) = ((1 - gamma2) * k2) / A3;
B(4,1) = ((1 - gamma1) * k1) / A4;

C = [1 0 0 0;
     0 1 0 0];
M = C;

Q1 = eye(2) * 1;
Q2 = eye(2) * 1;

% Solve CARE using icare
Q = M' * Q1 * M;   % Construct Q
R = Q2;            % R matrix
[S, ~, ~, ~] = icare(A, B, Q, R);

% Feedback gains
L = Q2\B'*S
Lr = inv(M*inv(B*L - A)*B)

%% Task 5: Performance analysis

%% Task 6: Robustness analysis

h10_i = 12.3;
h20_i = 12.8;
h30_i = 1.6;
h40_i = 1.4;
u_10_i = 3;
u_20_i = 3;
k1_i = 3.33;
k2_i = 3.35;

% G11
w = logspace(-4,0,50);      % frequency grid
l_I11 = zeros(size(w));    % max amplitude at each frequency

for i = 1:10000
    a1_i = a1 * (1 + 0.2 * (rand - 0.5));
    a2_i = a2 * (1 + 0.2 * (rand - 0.5));
    gamma1_i = 0.56 + (0.7 - 0.56) * rand;
    gamma2_i = 0.48 + (0.6 - 0.48) * rand;

    % Time constants
    T1_i = (A1/a1_i) * sqrt(2 * h10_i / g); 
    T2_i = (A2/a2_i) * sqrt(2 * h20_i / g);
    T3_i = (A3/a3) * sqrt(2 * h30_i / g);
    T4_i = (A4/a4) * sqrt(2 * h40_i / g); 

    c1_i = (T1_i * k1_i * kc) / A1;
    c2_i = (T2_i * k2_i * kc) / A2;

    G11_0_i = (gamma1_i * c1_i) / (T1_i*s + 1);
    deltaG = (G11_0_i - G11_M) / G11_M;

    [m, ~, ~] = bode(deltaG, w);
    m = squeeze(m);

    l_I11 = max(l_I11, m');   % update max at each frequency
end

figure;
loglog(w, l_I11, 'b', 'LineWidth', 2); hold on;
xlabel('\omega [rad/s]');
ylabel('|l_{I11}(i\omega)| (blue), |w_{I11}(i\omega)| (red)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
grid on;

w_I11 = 0.273 / (s/0.019 + 1) * (s/0.0256 + 1);
[m, ~, ~] = bode(w_I11, w);
m = squeeze(m);
loglog(w, m, 'r-', 'LineWidth', 2); hold off;

% G12
l_I12 = zeros(size(w));

for i = 1:10000
    a1_i = a1 * (1 + 0.2 * (rand - 0.5));
    a2_i = a2 * (1 + 0.2 * (rand - 0.5));
    gamma1_i = 0.56 + (0.7 - 0.56) * rand;
    gamma2_i = 0.48 + (0.6 - 0.48) * rand;

    % Calculate time constants T_i
    T1_i = (A1/a1_i) * sqrt(2 * h10_i / g); 
    T2_i = (A2/a2_i) * sqrt(2 * h20_i / g);
    T3_i = (A3/a3) * sqrt(2 * h30_i / g);
    T4_i = (A4/a4) * sqrt(2 * h40_i / g); 

    c1_i = (T1_i * k1_i * kc) / A1;
    c2_i = (T2_i * k2_i * kc) / A2;

    G12_0_i = ((1 - gamma2_i) * c1_i) / ((T1_i*s + 1) * (T3_i*s + 1)); 
    deltaG = (G12_0_i - G12_M) / G12_M;
    
    [m, ~, ~] = bode(deltaG, w);
    m = squeeze(m);

    l_I12 = max(l_I12, m');
end

figure;
loglog(w, l_I12, 'b', 'LineWidth', 2); hold on;
xlabel('\omega [rad/s]');
ylabel('|l_{I12}(i\omega)| (blue), |w_{I12}(i\omega)| (red)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
grid on;

w_I12 = 0.445 / (s/0.0172 + 1) * (s/0.0252 + 1);
[m, ~, ~] = bode(w_I12, w);
m = squeeze(m);
loglog(w, m, 'r-', 'LineWidth', 2); hold off;

% G21
l_I21 = zeros(size(w));

for i = 1:10000
    a1_i = a1 * (1 + 0.2 * (rand - 0.5));
    a2_i = a2 * (1 + 0.2 * (rand - 0.5));
    gamma1_i = 0.56 + (0.7 - 0.56) * rand;
    gamma2_i = 0.48 + (0.6 - 0.48) * rand;

    % Calculate time constants T_i
    T1_i = (A1/a1_i) * sqrt(2 * h10_i / g); 
    T2_i = (A2/a2_i) * sqrt(2 * h20_i / g);
    T3_i = (A3/a3) * sqrt(2 * h30_i / g);
    T4_i = (A4/a4) * sqrt(2 * h40_i / g); 

    c1_i = (T1_i * k1_i * kc) / A1;
    c2_i = (T2_i * k2_i * kc) / A2;

    G21_0_i = ((1 - gamma1_i) * c2_i) / ((T2_i*s + 1) * (T4_i*s + 1));
    deltaG = (G21_0_i - G21_M) / G21_M;
    [m, ~, ~] = bode(deltaG, w);
    m = squeeze(m);

    l_I21 = max(l_I21, m');
end

figure;
loglog(w, l_I21, 'b-', 'LineWidth', 2); hold on;
xlabel('\omega [rad/s]');
ylabel('|l_{I21}(i\omega)| (blue), |w_{I21}(i\omega)| (red)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
grid on;

w21 = 0.63 / (s/0.0135 + 1) * (s/0.018 + 1);
[m, ~, ~] = bode(w21, w);
m = squeeze(m);
loglog(w, m, 'r-', 'LineWidth', 2);

hold off;

% G22
l_I22 = zeros(size(w));

for i = 1:10000
    a1_i = a1 * (1 + 0.2 * (rand - 0.5));
    a2_i = a2 * (1 + 0.2 * (rand - 0.5));
    gamma1_i = 0.56 + (0.7 - 0.56) * rand;
    gamma2_i = 0.48 + (0.6 - 0.48) * rand;

    % Calculate time constants T_i
    T1_i = (A1/a1_i) * sqrt(2 * h10_i / g); 
    T2_i = (A2/a2_i) * sqrt(2 * h20_i / g);
    T3_i = (A3/a3) * sqrt(2 * h30_i / g);
    T4_i = (A4/a4) * sqrt(2 * h40_i / g); 

    c1_i = (T1_i * k1_i * kc) / A1;
    c2_i = (T2_i * k2_i * kc) / A2;

    G22_0_i = (gamma2_i * c2_i) / (T2_i*s + 1);
    deltaG = (G22_0_i - G22_M) / G22_M;
    [m, p, w] = bode(deltaG, w);
    m = squeeze(m);

    l_I22 = max(l_I22, m');
end

figure;
loglog(w, l_I22, 'b-', 'LineWidth', 2); hold on;
xlabel('\omega [rad/s]');
ylabel('|l_{I22}(i\omega)| (blue), |w_{I22}(i\omega)| (red)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
grid on;

w_I22 = 0.275 / (s/0.014 + 1) * (s/0.019 + 1);
[m, ~, ~] = bode(w_I22, w);
m = squeeze(m);
loglog(w, m, 'r-', 'LineWidth', 2); hold off;

% io(1) = linio('R7014E_QuadTankSim_LQG/Mux5', 1, 'input');
% io(2) = linio('R7014E_QuadTankSim_LQG/Sum4', 1, 'output');
% S = tf(linearize('R7014E_QuadTankSim_LQG', io));



% deltaG = [deltaG11, deltaG12;
%           deltaG21, deltaG22];
% 
% io(1) = linio('R7014E_QuadTankSim_LQG/Mux4', 1, 'input')  % y
% io(2) = linio('R7014E_QuadTankSim_LQG/Sum4', 1, 'output') % u
% F_y = linearize('R7014E_QuadTankSim_LQG', io);
% 
% io(1) = linio('R7014E_QuadTankSim_LQG/Sum4', 1, 'input')
% io(2) = linio('R7014E_QuadTankSim_LQG/Mux4', 1, 'output')
% G = linearize('R7014E_QuadTankSim_LQG', io);
% 
% L = G * F_y;
% T = feedback(L, eye(size(L)));
% T = d2c(T, 'tustin');
% 
% figure;
% sigma(deltaG * T);
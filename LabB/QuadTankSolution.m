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

% Tank water level initial values, tank 1-4
h1_init = 0;
h2_init = 0;
h3_init = 0;
h4_init = 0;

%% Operating Point Selection
h10_M = 12.3;
h20_M = 12.8;
h30_M = 1.6;
h40_M = 1.4;
u_10_M = 3;
u_20_M = 3;

h10_P = 12.6;
h20_P = 13.0;
h30_P = 4.8;
h40_P = 4.9;
u_10_P = 3.15;
u_20_P = 3.15;

k1_M = 3.33;
k2_M = 3.35;
gamma1_M = 0.7;
gamma2_M = 0.6;

k1_P = 3.14;
k2_P = 3.29;
gamma1_P = 0.43;
gamma2_P = 0.34;

gamma1 = gamma1_M;
gamma2 = gamma2_M;

k1 = k1_M;
k2 = k2_M;

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

%% Task 3
G11_num = cell2mat(G11_M.Numerator);
G11_den = cell2mat(G11_M.Denominator);
G11_inv = inv(G11_M);
lambda = 8;
Q11 = G11_inv/(lambda*s + 1)^2;
Q11_num = cell2mat(Q11.Numerator);
Q11_den = cell2mat(Q11.Denominator);

G22_num = cell2mat(G22_M.Numerator);
G22_den = cell2mat(G22_M.Denominator);
G22_inv = inv(G22_M);
lambda = 8;
Q22 = G22_inv/(lambda*s + 1)^2;
Q22_num = cell2mat(Q22.Numerator);
Q22_den = cell2mat(Q22.Denominator);

%% Task 4:

G_M_ss = ss(G_M);
A_size = size(G_M_ss.A);
B_size = size(G_M_ss.B);
C_size = size(G_M_ss.C);
Aa = [G_M_ss.A zeros(C_size(2), C_size(1));
     -G_M_ss.C zeros(C_size(1), C_size(1))];
Ba = [G_M_ss.B ; zeros(C_size(1), B_size(2))];
Ca = [G_M_ss.C, zeros(C_size(1), C_size(1))];

Q = eye(8);
R = 1;
Ka = lqr(ss(Aa,Ba,Ca,G_M_ss.D), Q, R);

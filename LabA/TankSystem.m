% (in cm)
At = 33;
Ao = 0.16;
g = 981;
Gain = 4; % Used for linearization 
T = 35;
h_0 = 3.2;
u0 = 1.07;

% Some constants
df_u_0_At = 7.1; % 7.54
Ts = 0.02;

% Equilibrium point
h_10 = 10; 
h_20 = 10;

% System equation
A = [(-Ao*g)/(At*sqrt(2*g*(h_10+h_0))) 0 ; 
    (Ao*g)/(At*sqrt(2*g*(h_10+h_0))) (-Ao*g)/(At*sqrt(2*g*(h_20+h_0)))];
B = [df_u_0_At/At; 0];
N = eye(2);
C = [0 1];
D = 0;

sys = ss(A,B,C,D);

% Discrete system
sys_d = c2d(sys,Ts);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

% Measurment Disturbance
Variance = 1e-6;

% Kalman Filter gain
Q = [Variance, 0; 
     0,    Variance];
R = 0.01;
[K, P] = dlqe(Ad, eye(size(Ad)), Cd, Q, R);

% Some mesurements
%u_measured = [0.4 0.8 1 1.1 1.2 1.3 1.9 2.5];
h_1_V = [0.54 0.965 1.34 1.72 2.01 2.43 2.86 3.215 3.62 4.01 4.39 4.76 5.17 5.51 5.91 6.28 6.66 7.05 7.49 7.83 8.2 8.5 9.3 9.6 9.95];
h_1 =[2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14.5 15];

h_2_V = [0.75 1.03 1.47 1.85 2.22 2.61 2.93 3.35 3.72 4.09 4.49 4.88 5.21 5.62 6.00 6.33 6.75  7.13 7.50 7.94 8.29 8.65 9.03 9.38 9.77];
h_2 =[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15];

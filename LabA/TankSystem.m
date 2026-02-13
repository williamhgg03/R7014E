% (in cm)
At = 33;
Ao = 0.16;
g = 981;
Gain = 4; % Used for linearization 
T = 35;
h_0 = 3.2;
Ts = 0.02;
Delay = 0; % in sec


% Equilibrium point
dq_0 = 7.28; % 7.54
u_0 = 1.39364;
h_10 = 8;
h_20 = 8;
Al_0 = 0;
x_hat_0 = [h_10; h_20];
x_hat_0_dist = [h_10; h_20; Al_0];

% System equation
sqrt1 = sqrt(2*g*(h_10+h_0));
sqrt2 = sqrt(2*g*(h_20+h_0));

A = [-Ao*g/At/sqrt1, 0;
     Ao*g/At/sqrt1, -Ao*g/At/sqrt2];
B = [dq_0/At; 0];
C = [0 1];
D = 0;

% Discrete system
Ad = eye(size(A)) + Ts * A;
Bd = Ts * B;
Cd = C;
Dd = D;

% Measurment Disturbance
Variance = 1e-5;

% Kalman Filter gain
Q = [Variance, 0; 
     0,    Variance];
R = 0.2;
[K, P] = dlqe(Ad, eye(size(Ad)), Cd, Q, R);

% Some mesurements
%u_measured = [0.4 0.8 1 1.1 1.2 1.3 1.9 2.5];

% Is the sensor fucked? measurment from upper tank does not seem to reflect
% reality...
h_1_V = [0.54 0.965 1.34 1.72 2.01 2.43 2.86 3.215 3.62 4.01 4.39 4.76 5.17 5.51 5.91 6.28 6.66 7.05 7.49 7.83 8.2 8.5 9.3 9.6 9.95];
h_1 =[2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14.5 15];

h_2_V = [0.75 1.03 1.47 1.85 2.22 2.61 2.93 3.35 3.72 4.09 4.49 4.88 5.21 5.62 6.00 6.33 6.75  7.13 7.50 7.94 8.29 8.65 9.03 9.38 9.77];
h_2 =[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15];

% Augmented system
Aa = [-(Ao+Al_0)*g/At/sqrt1, 0,              -sqrt1/At;
      Ao*g/At/sqrt1,         -Ao*g/At/sqrt2, 0;
      0,                     0,              0];
Ba = [dq_0/At; 0; 0];
Ca = [0 1 0];
Da = 0;

% Discrete system
Aad = eye(size(Aa)) + Ts * Aa;
Bad = Ts * Ba;
Cad = Ca;
Dad = Da;

% Kalman Filter gain with disturbance
Q_dist = [Variance 0 0; 
          0 Variance 0;
          0 0 1e-8];
[K_dist, P_dist] = dlqe(Aad, eye(size(Aad)), Cad, Q_dist, R);

% (in cm)
At = 33;
Ao = 0.16;
g = 981;
K = 4;
T = 35;
h_0 = 3.2;
df_u_0_At = 7.54;

% equilibrium point
h_10 = 8; 
h_20 = 8;

%System equation
A = [-Ao*g*(1/(At*sqrt(2*g*(h_10+h_0)))) 0 ; 
    Ao*g*(1/(At*sqrt(2*g*(h_10+h_0)))) Ao*g*(-1/(At*sqrt(2*g*(h_20+h_0))))];
B = [df_u_0_At/At; 0];
N = eye(2);
C = [0 1];
D = 0;

sys = ss(A,B,C,D);

u_measured = [0 1 1.9 2.5];
h_2_v_measured =[0 2.7 6.2 8.8];


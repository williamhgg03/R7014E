% (in cm)
At = 33;
Ao = 0.16;
g = 981;
K = 5;
T = 35;
h_0 = 3.2;

% equilibrium point
h_10 = 8; 
h_20 = 8;



A = [-Ao*g*(1/(At*sqrt(2*g*(h_10+h_0)))) 0 ; 
    Ao*g*(1/(At*sqrt(2*g*(h_10+h_0)))) -Ao*g*(-1/(At*sqrt(2*g*(h_20+h_0))))];
B = [6.95/36.2; 0];
N = eye(2);
C = [0 1];
D = 0;

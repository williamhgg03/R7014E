%% Script for AIRC control

%% Prepare the model
load airc

%% RGA analysis
% we need to remove the integrator
Gr1 = G(1,:);
Gr1n= minreal(Gr1*tf([1,0],[1]));
Grga= [Gr1n;G(2:3,:)];
Grga0=dcgain(Grga);
RGA = Grga0.*(inv(Grga0)');
fprintf(1,'-------------------------------------------------\n');
fprintf(1,'RGA analysis of the AIRC:\n');
display(RGA)
fprintf(1,'Most appropriate pairs:\n');
fprintf(1,'(y,u): (1,1), (2,3), (3,2)\n');


%%

G_Pm= (minreal(Grga));

P=gram(G_Pm,'c');
Q=gram(G_Pm,'o');

A=G_Pm.A;
B=G_Pm.B;
C=G_Pm.C;
D=G_Pm.D;

P_Test=0;

for i=1:3

    P_Test=P_Test+ gram(ss(A,B(:,i),C,D(:,i)),'c');
end
Q_Test=0;
for j=1:3

    Q_Test=Q_Test+ gram(ss(A,B ,C(j,:),D(j,:)),'o');
end



%%
Phi=0;
for i=1:3
    P_i=  gram(ss(A,B(:,i),C,D(:,i)),'c');
    for j=1:3
        Q_j= gram(ss(A,B ,C(j,:),D(j,:)),'o');
        phi(i,j)=trace(P_i*Q_j)/trace(P*Q);

    end
end

%% Problem 6.4
h = 1;
z = tf('z', 1);
Gz = (h^2)*(z+1)/(z-1)^2;
% Beta = 0
rlocus(Gz);
There are no values of K where CL system is valid with .
% Beta = 3/2
h = 1;
z = tf('z', 1);
Hz = (3/2)*(h/2)*(z-1)/(z+1);
G_new = feedback(Gz, Hz);
rlocus(G_new);
With , CL system is stable for .
K = 0.5;
NS = 4;
h = 1/NS;
beta = 3/2
z = tf('z', 1);
s = tf('s');
Gz = (h^2)*(z+1)/(z-1)^2;
Gs = 1/s^2;
Hz = beta*(h/2)*(z-1)/(z+1);
Hs = beta*s;
G_new = feedback(K*Gz, Hz);
Gs_new = feedback(K*Gs, Hs);
Tz = feedback(G_new, 1);
Ts = feedback(Gs_new, 1);
t = 0:h:10;
u = 1.5+heaviside(t);
figure;
lsim(Tz, u, t);
hold on;
plot(u, t);
legend('y',)

%% Nick's code %%
A = [0 1;0 0];
B = [0;1];
C = [1 0];
D = 0;
[n,d] = ss2tf(A,B,C,D);
Gs = tf(n,d); %continuous TF
Gz = c2d(Gs,1); %discretized TF
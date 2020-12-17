%% Problem 6.2 Rlocus %%
%% 5.1 %%
s = tf('s');
syms Ra La Kt Kv Bm;
Ra = 1;
La = .001;
Kt = 5;
Kv = 5;
Bm = 20;
J = 1;
A = [-Ra/La -Kv/La 0; Kt/J -Bm/J 0; 0 1 0];
B = [1/La; 0; 0];
C = [0 0 1];
[num, den] = ss2tf(A, B, C, 0);
A_ip = [0 0 1 0; 0 0 0 1; 0 11.5 -1.84e-3 -8.08e-5; 0 7.2 -8.16e-4 -9.18e-5];
B_ip = [0 0 .0184 .00816]';
C_ip = [1 0 0 0; 0 1 0 0];
C_ip_p = C_ip(1,:);
C_ip_theta = C_ip(2,:);
D_ip = [0];
sys_ip_p = ss(A_ip, B_ip, C_ip_p, D_ip);
sys_ip_theta = ss(A_ip,B_ip, C_ip_theta, D_ip);

Gs1 = 1/(8*s+1);
h1 = 1;
z1 = tf('z', h1);
Gz1 = c2d(Gs1, h1);
figure;
rlocus(Gz1);
grid;
Hz1_0 = (z1 - 0.8)/(z1 - 1);
figure;
rlocus(Gz1*Hz1_0);
grid;
K1 = 7;
Hz1 = K1*Hz1_0;
LG1 = Gz1*Hz1;
figure;
step(feedback(LG1, 1));
%% 5.2 %%
Gs2 = tf(num, den);
h2 = .005;
z2 = tf('z', h2);
Gz2 = c2d(Gs2, h2);
figure;
rlocus(Gz2);
grid;
Hz2_0 = (z2 - 0.8)/(z2 - .1);
figure;
rlocus(Gz2*Hz2_0);
grid;
K2 = 2000;
Hz2 = K2*Hz2_0;
LG2 = Gz2*Hz2;
figure;
step(feedback(LG2, 1));
%% 5.3 Maglev %%
Gs3_mag = -28.1/(s^2 - 1962); 
Gs3_mag_neg = -1 * Gs3_mag;
h3_mag = .004;
Gz3_mag_neg = c2d(Gs3_mag_neg, h3_mag);
Gz3_mag = -1*Gz3_mag_neg;
z3_mag = tf('z', h3_mag);
rlocus(Gz3_mag_neg);
grid;
Hz3_mag_0 = (z3_mag - 0.8)/(z3_mag + 1)
rlocus(Gz3_mag_neg*Hz3_mag_0);
grid;
K3_mag = -840;
Hz3_mag = K3_mag*Hz3_mag_0;
LG3_mag = Gz3_mag*Hz3_mag;
step(feedback(LG3_mag, 1));
%% 5.3 Inverted Pendulum position %%
close all;
Gs3_ip_p = tf(sys_ip_p);
h3_ipp = .07;
z3_ipp = tf('z', h3_ipp);
Gz3_ip_p = c2d(Gs3_ip_p, h3_ipp);
Gz3_ip_p_neg = -1*Gz3_ip_p;
%rlocus(Gz3_ip_p_neg);
figure;
rlocus(Gz3_ip_p);
grid;
%Hz3_ip_p_0 = (z3_ipp - 1.20655)*(z3_ipp-0.8)/(z3_ipp-1.1067)/(z3_ipp+0.2);
%Hz3_ip_p_0 = (z3_ipp-0.8)/(z3_ipp-1.2);
Hz3_ip_p_0 = (z3_ipp-0.8)/(z3_ipp+0.4)
figure;
rlocus(Gz3_ip_p*Hz3_ip_p_0);
grid;
K3_ipp = 10000;
%K3_ipp = 1;
Hz3_ip_p = K3_ipp*Hz3_ip_p_0;
LG3_ipp = Gz3_ip_p*Hz3_ip_p;
%pzmap(feedback(LG3_ipp, 1));
figure;
step(feedback(LG3_ipp, 1));
% Couldn't figure out
%% 5.3 Inverted Pendulum angle %%
close all;
Gs3_ip_theta = tf(sys_ip_theta);
h3_ipt = 0.0745;
z = tf('z', h3_ipt);
Gz3_ip_theta = c2d(Gs3_ip_theta, h3_ipt);
figure;
rlocus(Gz3_ip_theta);
grid;
Hz3_ipt_0 = (z-0.8)/(z-0.4);
figure;
rlocus(Gz3_ip_theta*Hz3_ipt_0);
grid;
K3_ipt = 4120;
LG3_ipt = K3_ipt*Hz3_ipt_0*Gz3_ip_theta;
figure;
step(feedback(LG3_ipt, 1));
% Hard to deal with unstable poles using rlocus method.
%% 5.4 %%
close all;
Gs4 = 1/(1+.08*s)^2;
h4 = 0.01;
z = tf('z', h4);
Gz4 = c2d(Gs4, h4);
figure;
rlocus(Gz4);
grid;
Hz4 = (z-0.8)/(z-1);
figure;
rlocus(Gz4*Hz4);
grid;
K4 = 0.9;
t = 0:h4:15;
v = 6*t;
LG4 = Gz4*Hz4*K4;
Tz = feedback(LG4, 1);
figure;
step(Tz)
figure;
lsim(Tz, v, t);
figure;
margin(LG4);
% CAN BE DONE BUT OVERSHOOT IS VERY HIGH
% Phase margin is about 30 degrees. 
%% 5.5 %%
close all;
Gs5 = tf([5],[0.05 1 0]);
h5 = 0.01;
z5 = tf('z', h5);
Gz5 = c2d(Gs5, h5);
figure;
rlocus(Gz5);
grid;
Hz5 = (z5 - 0.6)/(z5 + 0.2);
figure;
rlocus(Gz5*Hz5);
grid;
Kv = 1/h5*evalfr(Gz5*Hz5*(z5-1), 1);
K5 =64;
LG5 = Gz5*Hz5*K5;
Tz5 = feedback(LG5, 1);
figure;
margin(LG5);
t = 0: h5: 15;
u = sin(2*t);
figure;
lsim(Tz5, u, t);
%% 5.6 %%
close all;
Gs6 = tf([-1.117 3.1472], [1 4.6429 5.3821]);
h6 = 0.08;
z6 = tf('z', h6);
Gz6 = c2d(Gs6, h6);
figure;
rlocus(Gz6);
grid;
Hz6 = (z6-0.8)*(z6-0.7)/(z6-0.5)/(z6-1);
figure;
rlocus(Hz6*Gz6);
grid;
K6 = 1;
LG6 = K6*Gz6*Hz6;
Tz6 = feedback(LG6, 1);
figure;
step(Tz6);
% Settling time constraint was difficult to meet.
%% 5.7 %%
close all;
Gs7 = 1/s;
h7 = 0.1;
z7 = tf('z', h7);
Gz7 = c2d(Gs7, h7);
figure;
rlocus(Gz7);
grid;
Hz7 = (z7 - 0.4)/(z7-0.6)*(z7-0.7)/(z7-0.5);
figure;
rlocus(Gz7*Hz7);
grid;
K7 = 1;
LG7 = K7*Gz7*Hz7
figure;
margin(LG7);
% DIFFICULT USING RLOCUS
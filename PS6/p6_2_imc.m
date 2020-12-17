%% Problem 6.2 digital IMC %%
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
%% 5.1 %%
close all;
Gs1 = 1/(8*s+1);
h1 = 0.5;
z1 = tf('z', h1);
Gz1 = c2d(Gs1, h1);
Fz1 = (1-0.8)/(z1-0.8);
Qz1 = 1/Gz1*Fz1;
Hz1 = Qz1*1/(1-Gz1*Qz1);
Tz1 = feedback(Hz1*Gz1, 1);
step(Tz1);
title('Cruise Control Step response')
legend;
hold off;
%% 5.2 %%
close all;
Gs2 = tf(num, den);
h2 = .005;
z2 = tf('z', h2);
Gz2 = c2d(Gs2, h2);
b_minus = (1+1.659)*(1+0.0578); % zeros with negative real part
b_pluss = 1;
b_nm = 1;
Gz2_mod_inv = tf(poly([1 0.7937 0.007681]),b_minus*b_pluss*b_nm*4.044e-05, h2); 
Fz2 = (1-0.1)^3/(z2-0.1)^3;
Qz2 = Gz2_mod_inv*Fz2;
Hz2 = Qz2/(1-Gz2*Qz2);
figure;
margin(Hz2*Gz2);
%title("DC motor compensated stability margins");
figure;
step(feedback(Hz2*Gz2, 1));
title("DC motor compensated step response");
%% 5.3 Maglev %%
close all;
Gs3_mag = -28.1/(s^2 - 1962); 
Gs3_mag_neg = -1 * Gs3_mag;
h3_mag = .004;
Gz3_mag_neg = c2d(Gs3_mag_neg, h3_mag);
Gz3_mag = -1*Gz3_mag_neg;
z3_mag = tf('z', h3_mag);
b_minus = 2;
Gz3_mag_mod_inv = (z3_mag - 1.194)*(z3_mag - 0.8376)/(b_minus*-0.00022539);
Fz3_mag = (1-0.1)^2/(z3_mag-0.1)^2;
Qz3_mag = Gz3_mag_mod_inv*Fz3_mag;
Hz3_mag = Qz3_mag/(1-Gz3_mag*Qz3_mag);
LG3_mag = Gz3_mag*Hz3_mag;
figure;
margin(LG3_mag);
figure;
step(feedback(LG3_mag, 1));
title('Maglev compensated step response');
% Unstable system, will come back to it if time permits
%% 5.3 Inverted Pendulum position %%
close all;
Gs3_ip_p = tf(sys_ip_p);
h3_ipp = .07;
z3_ipp = tf('z', h3_ipp);
Gz3_ip_p = c2d(Gs3_ip_p, h3_ipp);
% Will come back to it if time permits
% RHP pole. Hard to do anything if we can't place poles where we want
%% 5.3 Inverted Pendulum angle %%
close all;
Gs3_ip_theta = tf(sys_ip_theta);
h3_ipt = 0.07;
z = tf('z', h3_ipt);
Gz3_ip_theta = c2d(Gs3_ip_theta, h3_ipt);
% Will come back if time permits
% Same problem
%% 5.4 %%
close all;
Gs4 = 1/(1+.08*s)^2;
h4 = 0.01;
z = tf('z', h4);
Gz4 = c2d(Gs4, h4);
b_minus = (1+0.92);
Gz4_mod_neg = (z-0.8825)^2/(b_minus*0.007191);
alpha = 0.5;
Fz4 = (1-alpha)^2/(z-alpha)^2;
Qz4 = Gz4_mod_neg*Fz4;
Hz4 = Qz4/(1-Gz4*Qz4);
t = 0:h4:15;
v = 6*t;
LG4 = Hz4*Gz4;
Tz = feedback(LG4, 1);
figure;
lsim(Tz, v, t);
figure;
margin(LG4);
%% 5.5 %%
close all;
Gs5 = tf([5],[0.05 1 0]);
h5 = 0.01;
z5 = tf('z', h5);
Gz5 = c2d(Gs5, h5);
b_minus = 1+0.9355;
Gz5_mod_neg = (z5-1)*(z5-0.8187)/(b_minus*0.0046827);
alpha5 = 0.05;
Fz5 = (1-alpha5)^2/(z5-alpha5)^2;
Qz5 = Gz5_mod_neg*Fz5;
Hz5 = Qz5/(1-Qz5*Gz5);
LG5 = Gz5*Hz5;
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
b_nm = -1.257*(z6-1/1.257);
Gz6_mod_neg = (z6-0.8361)*(z6-0.825)/(-0.065308*b_nm);
alpha6 = 0.5;
Fz6 = (1-alpha)/(z6-alpha);
Qz6 = Gz6_mod_neg*Fz6;
Hz6 = Qz6/(1-Qz6*Gz6);
LG6 = Hz6*Gz6;
Tz6 = feedback(LG6, 1);
figure;
margin(Hz6);
figure;
step(Tz6);
% Settling time constraint met 
%% 5.7 %%
%close all;
Gs7 = 1/s;
h7 = 0.03;
z7 = tf('z', h7);
Gz7 = c2d(Gs7, h7);
Gz7_mod_neg = (z7-1)/h7;
alpha7 = 0.9;
k = 4;
Fz7 = (1-alpha7)^k/(z7-alpha7)^k;
Qz7 = Gz7_mod_neg*Fz7;
Hz7 = Qz7/(1-Qz7*Gz7);
LG7 = Gz7*Hz7;
%figure;
margin(LG7);
%crossover frequency requirement difficult to meet
%required tuning of three things: k, alpha, and h.
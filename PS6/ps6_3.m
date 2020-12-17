%% Problem 6.3 %%
close all;
h = 1 % random h doesn't matter;
z = tf('z', h);
Ga = tf([1 0],[1 1.2 0 -0.2],h);
figure;
rlocus(Ga);
Gb = tf([-1 2.3], [1 -.1 -.6], h);
figure;
rlocus(Gb);

%% Nick's code %%
p1 = tf([0.96 0.6194 -2.469536 0.924768 0],[-1 -0.5208 0.864],-1);
roots([-1 -0.5208 0.864])
figure(1)
rlocus(p1);
p2=tf([22.6941 -22.8988 -0.2945 3.1044 0.384],[-5.29 2.76 0.64],-1);
roots([-5.29 2.76 0.64])
figure(2)
rlocus(p2);
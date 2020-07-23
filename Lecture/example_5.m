%% example 5-2
% zero pole gain
sys = tf([1 -2], [2 3 5]);
[z, p, k] = zpkdata(sys);

%
figure(1)
pzmap(sys)
grid on
%% example 5-3
close; clear; clc;
H1 = tf(1, [1 1]);
H2 = tf(10, [1 10]);
H3 = zpk([], [-10 -10], 100);
H4 = zpk([], [-1 -10], 10);

dt = 0.01;
t = 0:dt:8;
figure(2);
hold on;
step(H1, t);
step(H2, t);
step(H3, t);
step(H4, t);
hold off;

legend('H1', 'H2', 'H3', 'H4');
%% example 5-4
close; clear; clc;
H1 = tf(1, [1 1]);
H5 = zpk([], [-1 -1i 1i], 1);
H6 = zpk([], [-1 -1i -1i 1i 1i], 1);

dt = 0.01;
t = 0:dt:20;
figure(3);
hold on;
step(H1, t);
step(H5, t);
step(H6, t);
hold off;

legend('H1', 'H5', 'H6');
%% example 5-5
close; clear; clc;
H1 = zpk([], [-1 -2 -2 -0.5+4j -0.5-4j], 180);
H2 = zpk(-3, [-1 -2 -2 -0.5+4j -0.5-4j], 60);
H3 = zpk([-3 3], [-1 -2 -2 -0.5+4j -0.5-4j], -20);
H4 = zpk([3 -1-2j -1+2j], [-1 -2 -2 -0.5+4j -0.5-4j], -12);

dt = 0.01;
t = 0:dt:10;
figure(4);
subplot(2, 4, [1 2 5 6]);
hold on;
step(H1, t);
step(H2, t);
step(H3, t);
step(H4, t);
hold off;
grid on;
legend('H1', 'H2', 'H3', 'H4');

subplot(2, 4, 3);
pzmap(H1);
title('H1 pole-zero map');
% grid on;

subplot(2, 4, 4);
pzmap(H2);
title('H2 pole-zero map');
% grid on;

subplot(2, 4, 7);
pzmap(H3);
title('H3 pole-zero map');
% grid on;

subplot(2, 4, 8);
pzmap(H4);
title('H4 pole-zero map');
% grid on;
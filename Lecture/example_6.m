%% example 6-1
close; clear; clc;

R = 100;
C = 0.01;
sys = tf(1, [1 1]);

figure(1);
bode(sys);

%% example 6-2
dt = 0.01;

w1 = 1/R/C;
w2 = w1/10;
w3 = w1*10;
f1 = w1/2/pi;
f2 = w2/2/pi;
f3 = w3/2/pi;
t1 = 0:dt:(8/f1);
t2 = 0:dt:(8/f2);
t3 = 0:dt:(8/f3);

u1 = sin(w1*t1);
u2 = sin(w2*t2);
u3 = sin(w3*t3);

y1 = lsim(sys, u1, t1);
y2 = lsim(sys, u2, t2);
y3 = lsim(sys, u3, t3);

figure(2);
subplot(3, 1, 1);
plot(t1, u1, t1, y1);
legend('Input', 'Output');
title('sin(\omegat)');
xlabel('Time [sec]');
ylabel('Amplitude [voltage]');
axis tight;

subplot(3, 1, 2);
plot(t2, u2, t2, y2);
legend('Input', 'Output');
title('sin(\omega/10 t)');
xlabel('Time [sec]');
ylabel('Amplitude [voltage]');
axis tight;

subplot(3, 1, 3);
plot(t3, u3, t3, y3);
legend('Input', 'Output');
title('sin(10\omegat)');
xlabel('Time [sec]');
ylabel('Amplitude [voltage]');
axis tight;
%% Example 7-1
clear all

dt = 0.1;
t = 0:dt:100;

raw_data = zeros(size(t));
avg_data = zeros(size(t));

for n = 1:length(t)
    data = 10 + 2*randn(1);
    raw_data(n) = data;
%     raw_data = [raw_data data];
    avg_data(n) = AvgFilter(data);
end

figure(1)
plot(t, raw_data, 'k:', t, avg_data);

%% Example 7-2
clear all
load('IR_data.mat');

t = ir_data(:, 1);
raw_data = ir_data(:, 2);

avg_data = zeros(size(t));

for n = 1:length(t)
    data = raw_data(n);
    avg_data(n) = MovAvgFilter(data);
end

figure(2)
plot(t, raw_data, 'k:', t, avg_data);

%% Example 7-3
clear all
load('IR_data.mat');

t = ir_data(:, 1);
raw_data = ir_data(:, 2);

avg_data = zeros(size(t));

for n = 1:length(t)
    data = raw_data(n);
    avg_data1(n) = MovAvgFilter(data, 10);
end

for n = 1:length(t)
    data = raw_data(n);
    avg_data2(n) = MovAvgFilter(data, 2);
end

for n = 1:length(t)
    data = raw_data(n);
    avg_data3(n) = MovAvgFilter(data, 40);
end

figure(3)
plot(t, raw_data, 'k:');
hold on;
plot(t, avg_data1);
plot(t, avg_data2);
plot(t, avg_data3);
hold off;
legend('Raw Data', 'N=10', 'N=2', 'N=40');


%% Example 7-4
load('ir_data.mat');

t = ir_data(:, 1);
u = ir_data(:, 2);
wc = 4;
G0 = 1;

n = 1;
k = 1:n;
p = wc*exp(1j*(2*k+n-1)/2/n*pi);
bw_filt1 = zpk([], p, G0*wc^n);

n = 2;
k = 1:n;
p = wc*exp(1j*(2*k+n-1)/2/n*pi);
bw_filt2 = zpk([], p, G0*wc^n);

n = 3;
k = 1:n;
p = wc*exp(1j*(2*k+n-1)/2/n*pi);
bw_filt3 = zpk([], p, G0*wc^n);

n = 4;
k = 1:n;
p = wc*exp(1j*(2*k+n-1)/2/n*pi);
bw_filt4 = zpk([], p, G0*wc^n);

y1 = lsim(bw_filt1, u, t);
y2 = lsim(bw_filt2, u, t);
y3 = lsim(bw_filt3, u, t);
y4 = lsim(bw_filt4, u, t);

w = logspace(-2, 4, 1000);
figure(4);
bode(bw_filt1, w);
hold on;
bode(bw_filt2, w);
bode(bw_filt3, w);
bode(bw_filt4, w);
hold off;
grid on;

figure(5);
plot(t, u, 'k:', t, y1, t, y2, t, y3, t, y4);
grid on;
legend('Raw Signal', '1st Order', '2nd Order', '3rd Order', '4th Order');

%% Example 7-5
dt = 0.001;
t = 0:dt:20;
u = sin(t); % 1rad/sec sin(2*pi*f*t)

alpha = 0.01;
sys_integ = tf(1, [1 0]);
sys_diff = tf([1 0], [alpha 1]);

y_integ = lsim(sys_integ, u, t) - 1;
y_diff = lsim(sys_diff, u, t);

figure(6);
plot(t, u, 'k:', t, y_integ, t, y_diff);
legend('sine', 'integ', 'diff');

%% Example 7-6

dt = 0.1;
t = 0:dt:20;
u = sin(t);

% backward
y_bw = zeros(size(t));
for n = 2:length(t)
    y_bw(n) = (u(n) - u(n-1))/dt;
end

%forward
y_fw = zeros(size(t));
for n = 1:length(t)-1
    y_fw(n) = (u(n+1) - u(n))/dt;
end

% trapezoid
y_tp1 = zeros(size(t));
for n = 2:length(t)
    y_tp1(n) = 2/dt*(u(n)-u(n-1)) - y_tp1(n-1);
end

y_tp2 = zeros(size(t));
y_tp2(1) = 1;
for n = 2:length(t)
    y_tp2(n) = 2/dt*(u(n)-u(n-1)) - y_tp2(n-1);
end

figure(7);
subplot(1, 2, 1);
plot(t, u, ':');
hold on;
plot(t, y_bw);
plot(t, y_fw);
plot(t, y_tp2);
hold off;
legend('sin', 'backward', 'forward', 'trapz one init');

subplot(1, 2, 2);
plot(t, u, ':');
hold on;
plot(t, y_bw);
plot(t, y_fw);
plot(t, y_tp1);
hold off;
legend('sin', 'backward', 'forward', 'trapz zero init');
%% Example 7-7

dt = 0.1;
t = 0:dt:20;
u = sin(t);

y = -cos(t);
% backward
y_bw = zeros(size(t));
y_bw(1) = -1;
for n = 2:length(t)
    y_bw(n) = y_bw(n-1) + u(n-1)*dt;
end

%forward
y_fw = zeros(size(t));
y_fw(1) = -1;
for n = 2:length(t)
    y_fw(n) = y_fw(n-1) + u(n)*dt;
end

% trapezoid
y_tp = zeros(size(t));
y_tp(1) = -1;
for n = 2:length(t)
    y_tp(n) = y_tp(n-1) + dt/2*(u(n)+u(n-1));
end

figure(8);
plot(t, u, ':');
hold on;
plot(t, y);
plot(t, y_bw);
plot(t, y_fw);
plot(t, y_tp);
hold off;
legend('sin(t)', '-cos(t)', 'backward', 'forward', 'trapz');
%% Example 7-8
dt = 0.01;
t = 0:dt:10;
u = 1 < t & t < 2;

sys = tf(1, [1 1]);
y_lsim = lsim(sys, u, t);

y_z = zeros(size(t));
T = dt;
a0 = T/(2+T);
a1 = T/(2+T);
b1 = -(T-2)/(T+2);
for n = 2:length(t)
    y_z(n) = a0*u(n) + a1*u(n-1) + b1*y_z(n-1);
end

figure(9);
plot(t, u, ':', t, y_lsim, t, y_z);
legend('Input', 'Linear Simulation', 'Z-transform');
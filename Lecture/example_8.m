%% Example 8-1
dx = 0.5;
x = 0:dx:5;
y = 2*x+1 + randn(size(x));

Y = y';
X = [x' ones(size(x'))];
A_hat = (X'*X)\X'*Y;

figure(1);
plot(x, y, '*', x, 2*x+1, x, A_hat(1)*x+A_hat(2));
legend('Measured', 'Ideal', 'Estimated');

%% Example 8-2
dt = 0.01;
t = (0:dt:20)';
H = tf(1, [1 0.5 1.2]);
u = sin(t);
y = lsim(H, u, t) + 0.05*randn(size(t));

X = [sin(t) cos(t) ones(size(t))];
Y = y;
A_hat = X'*X\X'*Y;

Mag = sqrt(A_hat(1)^2 + A_hat(2)^2);
Phs = atan2(A_hat(2), A_hat(1));
Offset = A_hat(3);

figure(2);
plot(t, u, t, y, ':', t, Mag*sin(t+Phs)+Offset);
legend('Input', 'Measured', 'Estimated');

%% Example 8-3
% find steady-state response region
dt = 0.01;
t = (0:dt:40)';
u = ones(size(t));
H = tf(1, [1 0.5 1.2]);
y = lsim(H, u, t) + 0.05*randn(size(t));
tss = 15;

figure(3);
subplot(2, 1, 1);
plot(t, u, t, y);

% steady-state estimation
t = (0:dt:tss+20)';
u = sin(t);
y = lsim(H, u, t) + 0.05*randn(size(t));

X = [sin(t(t > tss)) cos(t(t > tss)) ones(size(t(t > tss)))];
Y = y(t > tss);
A_hat = X'*X\X'*Y;

Mag = sqrt(A_hat(1)^2 + A_hat(2)^2);
Phs = atan2(A_hat(2), A_hat(1));
Offset = A_hat(3);

subplot(2, 1, 2);
plot(t, u, t, y, ':', t, Mag*sin(t+Phs)+Offset);
legend('Input', 'Measured', 'Estimated');

%% Example 8-4
dt = 0.01;
t = (0:dt:40)';
w = logspace(-1, 1, 20);
H = tf(1, [1 0.5 1.2]);
tss = 15;

Mag = zeros(size(w));
Phs = zeros(size(w));
Offset = zeros(size(w));
for n = 1:length(w)
    u = sin(w(n)*t);
    y = lsim(H, u, t) + 0.05*randn(size(t));
    
    X = [sin(w(n)*t(t > tss)) cos(w(n)*t(t > tss)) ones(size(t(t > tss)))];
    Y = y(t > tss);
    A_hat = X'*X\X'*Y;
    
    Mag(n) = sqrt(A_hat(1)^2 + A_hat(2)^2);
    Phs(n) = atan2(A_hat(2), A_hat(1));
    Offset(n) = A_hat(3);
end
% Phs = unwrap(Phs);
[Mag_r, Phs_r] = bode(H, w);
Mag_r = Mag_r(:)';
Phs_r = Phs_r(:)';

figure(4);
subplot(2, 1, 1);
semilogx(w, 20*log10(Mag), '.', w, 20*log10(Mag_r)); % Magnitude
ylabel('Magnitude [dB]');
xlabel('Frequency [rad/sec]');
subplot(2, 1, 2);
semilogx(w, rad2deg(Phs), '.', w, Phs_r); % Phase
ylabel('Phase [\circ]');
xlabel('Frequency [rad/sec]');

[b, a] = invfreqs(Mag.*exp(1j*Phs), w, 0, 2)

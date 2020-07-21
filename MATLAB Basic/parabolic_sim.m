%% Motion Calculation
% init
theta = deg2rad(30);
v0 = 80*1000/3600; %[m/s]
h = 1.7; %[m]
g = 9.81; %[m/s^2]

% time scale
t_total = max(roots([g/2 -v0*sin(theta) -h]));
dt = 0.05;
t = 0:dt:t_total;
t = [t t_total];

% position calc
y = h + v0*sin(theta)*t - g/2*t.^2;
x = v0*cos(theta)*t;
%%
figure(1)
% data plot
subplot(2, 2, 3)
plot(t, x);
grid on;
title('x direction');
xlabel('Time [sec]');
ylabel('Distance [m]');
subplot(2, 2, 4);
plot(t, y);
grid on;
title('y direction');
xlabel('Time [sec]');
ylabel('Distance [m]');

% sim plot
subplot(2, 2, [1 2]);
n = 1;
ball = plot(x(n), y(n), 'ko');
hold on
grid on
ground = plot([-1 50], [0 0], 'k');
axis equal
xlim([-1 50]);
ylim([-1 10]);
xlabel('X Axis');
ylabel('Y Axis');

%%
tic
for n = 1:length(t)
    while toc < t(n)
    end
    title(['Simulation Time: ' sprintf('%.2f', t(n))])
    % update data of plot
    ball.XData = x(n);
    ball.YData = y(n);
    % update figure
    shg;
end
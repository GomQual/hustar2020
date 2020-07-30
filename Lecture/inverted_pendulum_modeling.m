%% Inverted Pendulum Cart Modeling
syms M m l b g J B x xp xpp th thp thpp F real
A = [m*l*cos(th) J+m*l^2; M+m m*l*cos(th)];
B = [m*g*l*sin(th)-B*thp; F-b*xp+m*l*sin(th)*thp^2];
A\B

%% Simulation Initial Value
M = 0.1;
m = 0.02;
l = 0.5;
b = 0.4;
g = 9.81;
J = m*l^2;
B = 0.001;

% initial value
dt = 0.001;
t = 0:dt:20;
xp = zeros(size(t));
thp = zeros(size(t));
x = zeros(size(t));
th = zeros(size(t));
F = zeros(size(t));
th(1) = deg2rad(10);

% Runge-Kutta Simulation
f = @(xp, thp, th, F) [[m*l*cos(th) J+m*l^2; M+m m*l*cos(th)]\[m*g*l*sin(th)-B*thp; F-b*xp+m*l*thp^2*sin(th)]; xp; thp];
for n = 1:length(t)-1
    k1 = f(xp(n), thp(n), th(n), F(n));
    k2 = f(xp(n)+dt*k1(1)/2, thp(n)+dt*k1(2)/2, th(n)+dt*k1(4)/2, F(n));
    k3 = f(xp(n)+dt*k2(1)/2, thp(n)+dt*k2(2)/2, th(n)+dt*k2(4)/2, F(n));
    k4 = f(xp(n)+dt*k3(1), thp(n)+dt*k3(2), th(n)+dt*k3(4), F(n));
    dy = dt/6*(k1+2*k2+2*k3+k4);
    xp(n+1) = xp(n) + dy(1);
    thp(n+1) = thp(n) + dy(2);
    x(n+1) = x(n) + dy(3);
    th(n+1) = th(n) + dy(4);
end

%% Simple Plot
figure(1);
subplot(2, 1, 1);
plot(t, rad2deg(th), t, rad2deg(thp));
title('Pole Motion');
xlabel('Time [sec]');
ylabel('Angle [\circ, \circ/sec]');
legend('\theta(t)', 'd\theta/dt');
subplot(2, 1, 2);
plot(t, x, t, xp);
title('Cart Motion');
xlabel('Time [sec]');
ylabel('Motion [m, m/sec]');
legend('x(t)', 'dx/dt');

%% Real-time Plot
figure(2);
subplot(4, 2, [1 2 3 4]);
n = 1;
circx = 0.02*cos(linspace(0, 2*pi, 10));
circy = 0.02*sin(linspace(0, 2*pi, 10));
ground = plot([-10 10], [0 0], 'k');
axes_handle = gca;
hold on
cart_size = 0.05;
cart = plot([x(n)-cart_size x(n)-cart_size x(n)+cart_size x(n)+cart_size x(n)-cart_size],...
    cart_size*[1 -1 -1 1 1], 'k');
rod = plot([x(n) x(n)+l*sin(th(n))], [0 l*cos(th(n))], 'k');
pen = plot(circx+x(n)+l*sin(th(n)), circy+l*cos(th(n)), 'k');
hold off
axis equal;
ylim([-0.6 0.6]);

time_range = 5;
subplot(4, 2, 5);
x_plot = plot(t(1:n), x(1:n));
axis auto;
title('x(t)');
xlabel('Time [sec]');
ylabel('Position [m]');
xlim([0 time_range]);
x_axes = gca;

subplot(4, 2, 6);
xp_plot = plot(t(1:n), xp(1:n));
axis auto;
title('x`(t)');
xlabel('Time [sec]');
ylabel('Velocity [m/sec]');
xlim([0 time_range]);
xp_axes = gca;

subplot(4, 2, 7);
th_plot = plot(t(1:n), rad2deg(th(1:n)));
axis auto;
title('\theta(t)');
xlabel('Time [sec]');
ylabel('Angle [\circ]');
xlim([0 time_range]);
th_axes = gca;

subplot(4, 2, 8);
thp_plot = plot(t(1:n), rad2deg(thp(1:n)));
axis auto;
title('\theta`(t)');
xlabel('Time [sec]');
ylabel('Angular Velocity [\circ/sec]');
xlim([0 time_range]);
thp_axes = gca;
%%
tic;
while n < length(t)
    while toc < t(n)
    end
    title(axes_handle, ['Simulation Time: ' sprintf('%.2f', t(n))])
    cart.XData = [x(n)-cart_size x(n)-cart_size x(n)+cart_size x(n)+cart_size x(n)-cart_size];
    rod.XData = [x(n) x(n)+l*sin(th(n))];
    rod.YData = [0 l*cos(th(n))];
    pen.XData = circx+x(n)+l*sin(th(n));
    pen.YData = circy+l*cos(th(n));
    
    x_plot.XData = t(1:n);
    x_plot.YData = x(1:n);
    xp_plot.XData = t(1:n);
    xp_plot.YData = xp(1:n);
    th_plot.XData = t(1:n);
    th_plot.YData = rad2deg(th(1:n));
    thp_plot.XData = t(1:n);
    thp_plot.YData = rad2deg(thp(1:n));
    if t(n) >= 5
        x_axes.XLim = [t(n)-5 t(n)];
        xp_axes.XLim = [t(n)-5 t(n)];
        th_axes.XLim = [t(n)-5 t(n)];
        thp_axes.XLim = [t(n)-5 t(n)];
    end
    figure(gcf);
    n = n + 100;
end 
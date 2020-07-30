%% Example 9-1

x = linspace(-2*pi, 2*pi, 10001);

ys = sin(x);
ys1 = sin(0) + cos(0)*x;
ys2 = sin(0) + cos(0)*x - sin(0)/prod(1:2)*x.^2;
ys3 = sin(0) + cos(0)*x - sin(0)/prod(1:2)*x.^2 - cos(0)/prod(1:3)*x.^3;
ys4 = sin(0) + cos(0)*x - sin(0)/prod(1:2)*x.^2 - cos(0)/prod(1:3)*x.^3 + sin(0)/prod(1:4)*x.^4;

yc = cos(x);
yc1 = cos(0) - sin(0)*x;
yc2 = cos(0) - sin(0)*x - cos(0)/prod(1:2)*x.^2;
yc3 = cos(0) - sin(0)*x - cos(0)/prod(1:2)*x.^2 + sin(0)/prod(1:3)*x.^3;
yc4 = cos(0) - sin(0)*x - cos(0)/prod(1:2)*x.^2 + sin(0)/prod(1:3)*x.^3 + cos(0)/prod(1:4)*x.^4;

figure(1);
subplot(2, 1, 1);
plot(x, ys, x, ys1, x, ys2, x, ys3, x, ys4);
xlim([-2*pi 2*pi]);
ylim([-1.1 1.1]);
legend('sin(x)', '1st', '2nd', '3rd', '4th');

subplot(2, 1, 2);
plot(x, yc, x, yc1, x, yc2, x, yc3, x, yc4);
xlim([-2*pi 2*pi]);
ylim([-1.1 1.1]);
legend('cos(x)', '1st', '2nd', '3rd', '4th');

%% Example 9-2

h = 0.01;
t = 0:h:5;

y_anal = exp(t);

y_df = zeros(size(t));
y_df(1) = 1;
for n = 2:length(t)
    y_df(n) = (2+h)/(2-h)*y_df(n-1);
end

y_euler = zeros(size(t));
y_euler(1) = 1;
for n = 2:length(t)
    y_euler(n) = (1+h)*y_euler(n-1);
end

figure(2);
plot(t, y_anal, ':', t, y_df, t, y_euler);
legend('Analytic', 'Digital Filter', 'Euler Method');

%% Example 9-3
% syms y(t)
% eqn = diff(y, t, 2) + diff(y, t) + y == 1;
% Dy = diff(y,t);
% cond = [y(0)==1, Dy(0)==1];
% ySol(t) = dsolve(eqn,cond)

h = 0.01;
t = 0:h:20;
y_anal = (2*3^(1/2)*exp(-t/2).*sin((3^(1/2)*t)/2))/3 + 1;

y_euler = zeros(2, length(t));
y_euler(:, 1) = [1; 1];
A = [-1 -1; 1 0];
for n = 2:length(t)
    y_euler(:, n) = y_euler(:, n-1) + h*(A*y_euler(:, n-1) + [1; 0]);
end

figure(3);
plot(t, y_anal, t, y_euler(2, :));
legend('Analytic', 'Euler Method');

%% Example 9-4
% init
h = 0.01;
t = 0:h:20;

% analytic
y_anal1 = exp(t);
y_anal2 = (2*3^(1/2)*exp(-t/2).*sin((3^(1/2)*t)/2))/3 + 1;

% euler method
y_euler1 = zeros(size(t));
y_euler1(1) = 1;
for n = 2:length(t)
    y_euler1(n) = y_euler1(n-1) + h*y_euler1(n-1);
end

y_euler2 = zeros(2, length(t));
y_euler2(:, 1) = [1; 1];
A = [-1 -1; 1 0];
B = [1; 0];
for n = 2:length(t)
    y_euler2(:, n) = y_euler2(:, n-1) + h*(A*y_euler2(:, n-1) + B);
end

% runge-kutta method
y_rk1 = zeros(size(t));
y_rk1(1) = 1;
for n = 2:length(t)
    k1 = y_rk1(n-1);
    k2 = y_rk1(n-1) + h*k1/2;
    k3 = y_rk1(n-1) + h*k2/2;
    k4 = y_rk1(n-1) + h*k3;
    y_rk1(n) = y_rk1(n-1) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

y_rk2 = zeros(2, length(t));
y_rk2(:, 1) = [1; 1];
A = [-1 -1; 1 0];
B = [1; 0];
for n = 2:length(t)
    k1 = A*y_rk2(:, n-1) + B;
    k2 = A*(y_rk2(:, n-1) + h*k1/2) + B;
    k3 = A*(y_rk2(:, n-1) + h*k2/2) + B;
    k4 = A*(y_rk2(:, n-1) + h*k3) + B;
    y_rk2(:, n) = y_rk2(:, n-1) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

figure(4);
subplot(2, 2, 1);
plot(t, y_anal1, ':', t, y_euler1, t, y_rk1);
subplot(2, 2, 2);
plot(t, y_anal2, ':', t, y_euler2(2, :), t, y_rk2(2, :));

subplot(2, 2, 3);
plot(t, y_euler1-y_anal1, t, y_rk1-y_anal1);
subplot(2, 2, 4);
plot(t, y_euler2(2, :)-y_anal2, t, y_rk2(2, :)-y_anal2);

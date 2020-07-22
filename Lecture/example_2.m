%% example 2-1
% time series
dt = 0.01;
t = 0:dt:10;
% analytic solution
y_anal = 2*exp(-t/2)+3;

% plot analytic solution
figure(1);
plot(t, y_anal);
hold on

% analytic solution by symbolic math
% using dsolve()
syms y(t)
eqn = 2*diff(y, t) + y == 3;
cond = y(0) == 5;
y_solv = dsolve(eqn, cond);
% evaluation of symbolic solution
dt = 0.01;
t = 0:dt:10;
y_syms = eval(y_solv);

plot(t, y_syms)
hold off
%% example 2-2
% analytic solution by symbolic math
% using ilaplace()
syms s
eqn = 1/(s+4);
y_syms = ilaplace(eqn);
% time series
dt = 0.01;
t = 0:dt:10;
% evaluation of symbolic solution
y_solv = eval(y_syms);
figure(2);
plot(t, y_solv);
%% example 4-1 (overdamped)

% analytic solution
a2 = 1;
a1 = 4;
a0 = 2;
A = 3;
x0 = 1;
Dx0 = 0;

sig1 = (-a1 + sqrt(a1^2-4*a2*a0))/(2*a2);
sig2 = (-a1 - sqrt(a1^2-4*a2*a0))/(2*a2);
B1 = (Dx0-(x0-A/a0)*sig2)/(sig1-sig2);
B2 = -(Dx0-(x0-A/a0)*sig1)/(sig1-sig2);

dt = 0.01;
t = 0:dt:10;
x_anal = B1*exp(sig1*t) + B2*exp(sig2*t) + A/a0;

figure(1);
subplot(2, 1, 1);
plot(t, x_anal);

% using symbolic functions
syms x(t)

Dx = diff(x, t);
D2x = diff(x, t, 2);
eqn = a2*D2x + a1*Dx + a0*x == A;
cond = [ x(0) == x0, Dx(0) == Dx0 ];
x_syms = dsolve(eqn, cond);

t = 0:dt:10;
x_eval = eval(x_syms);

subplot(2, 1, 2);
plot(t, x_eval);

%% example 4-2 (critically damped)

% analytic solution
a2 = 1;
a1 = 4;
a0 = 4;
A = 3;
x0 = 1;
Dx0 = 0;

sig = -a1/(2*a2);
D1 = x0 - A/a0;
D2 = Dx0-(D1)*sig;

dt = 0.01;
t = 0:dt:10;
x_anal = (D1 + D2*t).*exp(sig*t) + A/a0;

figure(2);
subplot(2, 1, 1);
plot(t, x_anal);

% using symbolic functions
syms x(t)

Dx = diff(x, t);
D2x = diff(x, t, 2);
eqn = a2*D2x + a1*Dx + a0*x == A;
cond = [ x(0) == x0, Dx(0) == Dx0 ];
x_syms = dsolve(eqn, cond);

t = 0:dt:10;
x_eval = eval(x_syms);

subplot(2, 1, 2);
plot(t, x_eval);

%% example 4-3 (underdamped)

% analytic solution
a2 = 1;
a1 = 4;
a0 = 10;
A = 3;
x0 = 1;
Dx0 = 0;

sig = -a1/(2*a2);
w = sqrt(4*a2*a0-a1^2)/(2*a2);
C1 = x0 - A/a0;
C2 = (Dx0-(C1)*sig)/w;

dt = 0.01;
t = 0:dt:10;
x_anal = exp(sig*t).*(C1*cos(w*t)+C2*sin(w*t)) + A/a0;

figure(3);
subplot(2, 1, 1);
plot(t, x_anal);

% using symbolic functions
syms x(t)

Dx = diff(x, t);
D2x = diff(x, t, 2);
eqn = a2*D2x + a1*Dx + a0*x == A;
cond = [ x(0) == x0, Dx(0) == Dx0 ];
x_syms = dsolve(eqn, cond);

t = 0:dt:10;
x_eval = eval(x_syms);

subplot(2, 1, 2);
plot(t, x_eval);
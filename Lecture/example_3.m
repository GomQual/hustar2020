%% example 3-1
% time series
dt = 0.01;
t = 0:dt:3;
% signal description
u = t < 2;
r = t .* (t <= 1);

% plot raw signals
figure(1);
subplot(2, 4, [1 2 5 6]);
plot(t, u, t, r);
legend('u(t)', 'r(t)');

% convolutions
y1 = conv(u, u);
y2 = conv(u, r);
y3 = conv(r, u);
y4 = conv(r, r);

% extended time series (to match data length)
t2 = (0:(length(y1)-1))*dt;

% plot each results
subplot(2, 4, 3);
plot(t2, y1);
title('u*u');
subplot(2, 4, 4);
plot(t2, y2);
title('u*r');
subplot(2, 4, 7);
plot(t2, y3);
title('r*u');
subplot(2, 4, 8);
plot(t2, y4);
title('r*r');

%% example 3-2
%%% convolution of impulse response
% time series
dt = 0.01;
t = 0:dt:10;

% impulse response h(t)
h = 1/2*exp(-0.5*t);
% step input u(t)
u = ones(size(t));
% convolution
y1 = conv(h, u)*dt;
% cut data length to match time series
y1 = y1(1:length(t));

% plot convolution
figure(2);
subplot(3, 1, 1);
plot(t, y1);
title('Convolution with Impulse Response');

%%% inverse laplace - transfer function
% symbolic inverse laplace
syms s
eqn = 1/(2*s+1)*(1/s);
y_syms = ilaplace(eqn);

% time series
dt = 0.01;
t = 0:dt:10;

% evaluation
y2 = eval(y_syms);
subplot(3, 1, 2);
plot(t, y2);
title('Transfer Function');

%%% using control systme toolbox
sys = tf(1, [2 1]); 
y3 = step(sys, t);
subplot(3, 1, 3);
plot(t, y3);
title('Using step() Function');
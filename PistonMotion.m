clc
clear all

bF = 0.09;
a = -0.1;
w0 = 4.87;
epsilon = 0.2;
omega = 0.9;
psi = 0;

b0 = 1 - bF;
bn1 = (1 - b0 + a)/2;
bp1 = (1 - b0 - a)/2;
a0 = epsilon / 2.09439;
T0 = 1 / (w0 * epsilon);
dt = 0.00001;

x = 0;
X = [];
X(1) = x;
T2 = [];
t=0;
T2(1) = t;

while t<40
    if t < 0.32254
        x = -3.254 * t^2 + 1.052*t - 0.0007626;
        t = t + dt;
        T2(end+1) = t;
        X(end+1) = x;
    else
        a1 = sqrt(b0) * cos(w0 * t);
        T = t / T0;
        a2 = sqrt(bn1) * cos(omega * T + psi) * cos(w0 * t) - sin(omega*T + psi) * sin(w0 * t);
        a3 = abs(sqrt(bp1)) * cos(-omega * T + psi) * cos(w0 * t) - sin(-omega*T + psi) * sin(w0 * t);

        x = a0 * (a1 + a2 + a3);
        t = t + dt;
        T2(end+1) = t;
        X(end+1) = x;
    end
end
size(X)
size(T2)
figure
plot(T2,X)
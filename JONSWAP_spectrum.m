clc
clear all

h = 2.0;
Tp = 2.0;
Hs  = 0.6;
gamma = 3.3;
alpha = 0.0624 / (0.23 + 0.0336 * gamma - (0.185/(1.9+gamma)));
% 
fp = (1/Tp);
fStart = 0.2;
fStop = 2.5;

sigma = 0.07;
beta = exp(-(fStart-fp)^2/(2*sigma^2*fp^2));
SfStart = alpha * Hs^2 * fp^4 * fStart^(-5) * gamma^beta * exp((-5/4)*(fp/fStart)^4)

beta = exp(-(fp-fp)^2/(2*sigma^2*fp^2));
Sfp = alpha * Hs^2 * fp^4 * fp^(-5) * gamma^beta * exp((-5/4)*(fp/fp)^4)

beta = exp(-(fStop-fp)^2/(2*sigma^2*fp^2));
SfEnd = alpha * Hs^2 * fp^4 * fStop^(-5) * gamma^beta * exp((-5/4)*(fp/fStop)^4)

N = 1500;
Sf = zeros(1,N);
F = zeros(1,N);
F2 = zeros(N,1);
for i=1:N
    f = fStart + i*((fStop - fStart)/N);
    F(i) = f;
    F2(i) = f / fp;
    if f < fp
        sigma = 0.07;
        beta = exp(-(f-fp)^2/(2*sigma^2*fp^2));
        Sf(i) = alpha * Hs^2 * fp^4 * f^(-5) * gamma^beta * exp((-5/4)*(fp/f)^4);
    else
        sigma = 0.09;
        beta = exp(-(f-fp)^2/(2*sigma^2*fp^2));
        Sf(i) = alpha * Hs^2 * fp^4 * f^(-5) * gamma^beta * exp((-5/4)*(fp/f)^4);
    end
end

Tp2 = 4/3;
Hs2  = 0.6;

fp2 = (1/Tp2);
fStart2 = 0.2;
fStop2 = 2.5;

sigma = 0.07;
beta2 = exp(-(fStart2-fp2)^2/(2*sigma^2*fp2^2));
SfStart2 = alpha * Hs2^2 * fp2^4 * fStart2^(-5) * gamma^beta2 * exp((-5/4)*(fp2/fStart2)^4)

beta2 = exp(-(fp2-fp2)^2/(2*sigma^2*fp2^2));
Sfp2 = alpha * Hs2^2 * fp2^4 * fp2^(-5) * gamma^beta2 * exp((-5/4)*(fp2/fp2)^4)

beta2 = exp(-(fStop2-fp2)^2/(2*sigma^2*fp2^2));
SfEnd2 = alpha * Hs2^2 * fp2^4 * fStop2^(-5) * gamma^beta2 * exp((-5/4)*(fp2/fStop2)^4)

N = 1500;
Sf2 = zeros(1,N);
F2 = zeros(1,N);
F22 = zeros(N,1);
for i=1:N
    f2 = fStart2 + i*((fStop2 - fStart2)/N);
    F2(i) = f2;
    F22(i) = f2 / fp2;
    if f2 < fp2
        sigma = 0.07;
        beta2 = exp(-(f2-fp2)^2/(2*sigma^2*fp2^2));
        Sf2(i) = alpha * Hs2^2 * fp2^4 * f2^(-5) * gamma^beta2 * exp((-5/4)*(fp2/f2)^4);
    else
        sigma = 0.09;
        beta2 = exp(-(f2-fp2)^2/(2*sigma^2*fp2^2));
        Sf2(i) = alpha * Hs2^2 * fp2^4 * f2^(-5) * gamma^beta2 * exp((-5/4)*(fp2/f2)^4);
    end
end
SfTotal = Sf + Sf2;
size(SfTotal)
plot(F2,SfTotal)
xlim([F2(1) F2(end)]);
grid on;
grid minor;
xlabel('$$f$$','interpreter','latex')
ylabel('$$S_f(m^2/s)$$','interpreter','latex')
title('JONSWAP theoretical spectrum','interpreter','latex')
f = gcf;
print(f,'-dpng','Spectrum1_3.png')
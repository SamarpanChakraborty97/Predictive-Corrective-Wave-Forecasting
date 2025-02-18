clc
clear all

fl = 0.8944;
fh = 4.4721;
fp = 2.4597;
Tp = 1/fp;
H_13 = 1;
g = 3.3;
a = abs((0.06238/(0.0230 + 0.0336*g - 0.185*(1.9+g)))*(1.094 - 0.01915 * log(g)));
Nf = 200;
delf = (fh - fl)/Nf;
Ab = 0.0039;
F = zeros(Nf,1);
Sf = zeros(Nf,1);
tb = 26.8328;
xb = 4.4;

F(1) = fl;
Sf(1) = a * H_13^2*Tp^(-4)*(F(1))^(-5)*exp(-1.25*(Tp*F(1))^-4)*g^(exp(-(F(1)/fp)-1)^2 / 2*0.07^2);

for i=2:Nf
    F(i) = F(i-1) + delf;
    if F(i)>=fp
        s = 0.09;
        Sf(i) = a * (H_13^2)*(Tp^(-4)*(F(i))^(-5))*exp(-1.25*(Tp*F(i))^-4)*(g^(exp(-(F(i)/fp)-1)^2 / 2*s^2));
    else
        s = 0.07;
        Sf(i) = a * (H_13^2)*(Tp^(-4)*(F(i))^(-5))*exp(-1.25*(Tp*F(i))^-4)*(g^(exp(-(F(i)/fp)-1)^2 / 2*s^2));
    end
end

a =  Ab * Sf/ sum(Sf);
K = (2*pi*F).^2/9.81;

t=0:1e-5:40;
Eta = zeros(length(t),1);
for i = 1:length(t)
    sum = 0;
    time = t(i);
    for j=1:Nf
        sum = sum + a(j)*cos(K(j)*(-xb) - 2*pi*F(j)*(time-tb));
    end
    Eta(i) = sum;
end
plot(t,Eta)


 
% for i=1:Nf
%     syms ki
%     eqn = (2*pi*F(i))^2 == 9.81*ki*tanh(ki*1.2);
%     vpasolve(eqn,ki)
%     K(i) = k;
% end
% K
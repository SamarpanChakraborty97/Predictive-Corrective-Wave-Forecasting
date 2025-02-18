clc
clear all;

a = readmatrix("Wave_Spectrum_076_march23_23_30.txt");

Nf = length(a(:,1));
Freqs = a(:,1);
Energy = a(:,3);
del_f = a(:,2);

figure;
f = gcf;
bar(Freqs, Energy, 1.0, 'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
grid on;
grid minor;
xlabel('Frequency (Hz)','interpreter','latex')
ylabel('Energy density ($m^{2}/Hz$)','interpreter','latex')
title('Spectrum at buoy location at given time','interpreter','latex')
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 15 7];
name2 = sprintf('Wavemaker initial spectrum_2_23_30.png');
print(f,name2,'-dpng','-r300');


Hs = 2.26;
d = 27.4;
Tp = 9.88;
fp = 1/Tp;
syms kp;
eqn = (2*pi*fp)^2 == 9.81*kp*tanh(kp*d);
kp = double(vpasolve(eqn));
lambda_p = 2*pi / kp;

Omegas = 2*pi.*Freqs;
K = zeros(Nf,1);
amps = zeros(Nf,1);
tran_funcs = zeros(Nf,1);

for i=1:Nf
    syms ki;
    eqn = (2*pi*Freqs(i))^2 == 9.81*ki*tanh(ki*d);
    k = double(vpasolve(eqn));
    K(i) = k;
    amps(i) = sqrt(2*Energy(i)*del_f(i));
    tran_funcs(i) = 2*(cosh(2*k*d)-1)/(2*k*d + sinh(2*k*d));
end

P(:,1) = Omegas;
P(:,2) = amps;
P(:,3) = tran_funcs;
P(:,4) = K';

writematrix(P,'InitP.txt')

% x_wavemaker_1 = 1.02 * lambda_p;
% x_wavemaker_2 = 0;
% x_wavemaker_3 = -3;
% x_wavemaker_4 = -20.0;

% x_buoy_1 = 1.02 * lambda_p;
% x_buoy_2 = 0;
% x_buoy_3 = 3;
% x_buoy_4 = 20.0;
% 
% %x_wavemaker = 1.02 * lambda_p;
% dt = 0.01;
% t = 0;
% 
% disp_1(1) = 0;
% disp_2(1) = 0;
% disp_3(1) = 0;
% disp_4(1) = 0;
% T(1) = 0;
% 
% i=2;
% while t<=300
% %     displacement1 = 0;
% %     displacement2 = 0;
% %     displacement3 = 0;
% %     displacement4 = 0;
%     
%     vel1 = 0;
%     vel2 = 0;
%     vel3 = 0;
%     vel4 = 0;
%     
%     for j=1:Nf
% %         displacement1 = displacement1 + dt * (Omegas(j)/tran_funcs(j))*amps(j) * cos(K(j)*x_wavemaker_1 - Omegas(j)*t);
% %         displacement2 = displacement2 + dt * (Omegas(j)/tran_funcs(j))*amps(j) * cos(K(j)*x_wavemaker_2 - Omegas(j)*t);
% %         displacement3 = displacement3 + dt * (Omegas(j)/tran_funcs(j))*amps(j) * cos(K(j)*x_wavemaker_3 - Omegas(j)*t);
% %         displacement4 = displacement4 + dt * (Omegas(j)/tran_funcs(j))*amps(j) * cos(K(j)*x_wavemaker_4 - Omegas(j)*t);
%         
% %         vel1 = vel1 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(K(j)*(x_wavemaker_1) - Omegas(j)*t);
% %         vel2 = vel2 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(K(j)*x_wavemaker_2 - Omegas(j)*t);
% %         vel3 = vel3 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(K(j)*x_wavemaker_3 - Omegas(j)*t);
% %         vel4 = vel4 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(K(j)*x_wavemaker_4 - Omegas(j)*t);
%         
%         vel1 = vel1 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(Omegas(j)*t);
%         vel2 = vel2 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(Omegas(j)*t);
%         vel3 = vel3 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(Omegas(j)*t);
%         vel4 = vel4 + (Omegas(j)/tran_funcs(j))* amps(j) * cos(Omegas(j)*t);
%     end
% %     disp_1(i) = displacement1;
% %     disp_2(i) = displacement2;
% %     disp_3(i) = displacement3;
% %     disp_4(i) = displacement4;
%     
%     disp_1(i) = disp_1(i-1) + dt * vel1;
%     disp_2(i) = disp_2(i-1) + dt * vel2;
%     disp_3(i) = disp_3(i-1) + dt * vel3;
%     disp_4(i) = disp_4(i-1) + dt * vel4;
%     
%     T(i) = t;
%     t = t+dt;
%     i = i+1;
% end

% figure;
% f = gcf;
% 
% subplot(2,2,1)
% plot(T, disp_1, 'b-','Linewidth',0.5);
% grid on;
% grid minor;
% xlim([0 300]);
% xlabel('time (in seconds)','interpreter','latex')
% ylabel('$$\eta (m)$$','interpreter','latex')
% title('x=-170.46','interpreter','latex')
% 
% subplot(2,2,2)
% plot(T, disp_2, 'k-','Linewidth',0.5);
% grid on;
% grid minor;
% xlim([0 300]);
% xlabel('time (in seconds)','interpreter','latex')
% ylabel('$$\eta (m)$$','interpreter','latex')
% title('x=0','interpreter','latex')
% 
% subplot(2,2,3)
% plot(T, disp_3, 'r-','Linewidth',0.5);
% grid on;
% grid minor;
% xlim([0 300]);
% xlabel('time (in seconds)','interpreter','latex')
% ylabel('$$\eta (m)$$','interpreter','latex')
% title('x=-3.0','interpreter','latex')
% 
% subplot(2,2,4)
% plot(T, disp_4, 'm-','Linewidth',0.5);
% grid on;
% grid minor;
% xlim([0 300]);
% xlabel('time (in seconds)','interpreter','latex')
% ylabel('$$\eta (m)$$','interpreter','latex')
% title('x=-20.0','interpreter','latex')
% 
% f.PaperUnits = 'inches';
% f.PaperPosition = [0 0 15 7];
% name2 = sprintf('Wavemaker displacement.png');
% print(f,name2,'-dpng','-r300');
% 
% close all;
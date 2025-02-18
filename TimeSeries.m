%%%%%%%%%%%%%CREATING THE DISCRETIZED FREQUENCY SPECTRUM%%%%%%%%%%%%%%%
if 1    
    clc
    clear all

    h = 2.0;
    Tp = 2.0;
    Hs = 0.6;
    gamma = 3.3;
    alpha = 0.0624 / (0.23 + 0.0336 * gamma - (0.185/(1.9+gamma)));

    fp = (1/Tp);
    fStart = 0.2;
    fStop = 2.0;

    N = 128;
    delF = (fStop - fStart)/ N;
    f = zeros(N,1);
    f2 = f;
    k = f; L = f; S = f; w = f;B = f;H = f;
    for i = 1:N
        f(i) = fStart + i*delF - delF/2;
        f2(i) = f(i)/fp;
        w(i) = 2* pi * f(i);
        k(i) = w(i) * w(i)/9.81;
        L(i) = 2*pi/ k(i);
        B(i) = 2 * sinh(k(i)*h) * sinh(k(i)*h)/((sinh(k(i)*h) * cosh(k(i)*h))+k(i)*h);
        if f(i) < fp
            sigma = 0.07;
            beta = exp(-(f(i)-fp)^2/(2*sigma^2*fp^2));
            S(i) = alpha * Hs^2 * fp^4 * f(i)^(-5) * gamma^beta * exp((-5/4)*(fp/f(i))^4);
        else
            sigma = 0.09;
            beta = exp(-(f(i)-fp)^2/(2*sigma^2*fp^2));
            S(i) = alpha * Hs^2 * fp^4 * f(i)^(-5) * gamma^beta * exp((-5/4)*(fp/f(i))^4);
        end
        H(i) = 2 * sqrt(2 * S(i) * delF);
    end
    del = (2*pi).*rand(N,1);
    
    Matrix = zeros(N,4);
    Matrix(:,1) = H;
    Matrix(:,2) = B;
    Matrix(:,3) = w;
    Matrix(:,4) = del;
    
    bar(f2,S,1,'FaceColor',[0.5 0.6 0.8])
    xlim([f2(1) f2(end)]);
    grid on;
    grid minor;
    xlabel('$$f/f_p$$','interpreter','latex')
    ylabel('$$S_f(m^2/s)$$','interpreter','latex')
    title('JONSWAP experimental spectrum','interpreter','latex')
    f = gcf;
    print(f,'-dpng','Spectrum2.png')
end
writematrix(Matrix,'ParameterValues.txt','Delimiter',' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if 1    
    clc
    clear all

    h = 2.0;
    Tp2 = 4/3;
    Hs = 0.6;
    gamma = 3.3;
    alpha = 0.0624 / (0.23 + 0.0336 * gamma - (0.185/(1.9+gamma)));

    fp2 = (1/Tp2);
    fStart2 = 0.25;
    fStop2 = 2.5;

    N = 128;
    delF2 = (fStop2 - fStart2)/ N;
    f2 = zeros(N,1);
    f22 = f2;
    k2 = f2; L2 = f2; S2 = f2; w2 = f2;B2 = f2;H2 = f2;
    for i = 1:N
        f2(i) = fStart2 + i*delF2 - delF2/2;
        f22(i) = f2(i)/fp2;
        w2(i) = 2* pi * f2(i);
        k2(i) = w2(i) * w2(i)/9.81;
        L2(i) = 2*pi/ k2(i);
        B2(i) = 2 * sinh(k2(i)*h) * sinh(k2(i)*h)/((sinh(k2(i)*h) * cosh(k2(i)*h))+k2(i)*h);
        if f2(i) < fp2
            sigma = 0.07;
            beta2 = exp(-(f2(i)-fp2)^2/(2*sigma^2*fp2^2));
            S2(i) = alpha * Hs^2 * fp2^4 * f2(i)^(-5) * gamma^beta2 * exp((-5/4)*(fp2/f2(i))^4);
        else
            sigma = 0.09;
            beta2 = exp(-(f2(i)-fp2)^2/(2*sigma^2*fp2^2));
            S2(i) = alpha * Hs^2 * fp2^4 * f2(i)^(-5) * gamma^beta2 * exp((-5/4)*(fp2/f2(i))^4);
        end
        H2(i) = 2 * sqrt(2 * S2(i) * delF2);
    end
    del2 = (2*pi).*rand(N,1);
    
    Matrix2 = zeros(N,4);
    Matrix2(:,1) = H2;
    Matrix2(:,2) = B2;
    Matrix2(:,3) = w2;
    Matrix2(:,4) = del2;
    
    bar(f22,S2,1,'FaceColor',[0.5 0.6 0.8])
    xlim([f22(1) f22(end)]);
    grid on;
    grid minor;
    xlabel('$$f/f_p$$','interpreter','latex')
    ylabel('$$S_f(m^2/s)$$','interpreter','latex')
    title('JONSWAP experimental spectrum','interpreter','latex')
    f = gcf;
    print(f,'-dpng','Spectrum2_2.png')
end
writematrix(Matrix2,'ParameterValues2.txt','Delimiter',' ')
%%%%%%%%%%%%%%%CREATING THE TIME SERIES DATA%%%%%%%%%%%%%
if  1
t = 0;
dt = 1 * 10^(-5);
T = [];
P = [];
while t<30
    sum = 0;
    for i = 1:N
        sum = sum + H(i)/(2*B(i)) * cos(w(i)*t + del(i)+pi/2);
    end
    T(end+1) = t;
    P(end+1) = sum;
    t = t + dt;
end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
plot(T,P,'r','Linewidth',2)
xlim([T(1) T(end)]);
grid on;
grid minor;
xlabel('$$t(secs)$$','interpreter','latex')
ylabel('$$S(t) (metres)$$','interpreter','latex')
title('Piston movement time series','interpreter','latex')
f = gcf;
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 10 3];
set(gca,'FontName','Times','FontSize',12);
print(f,'-dpng','Piston movement time series.png')
saveas(f,'Piston movement time series2.png')

end
%%%%%%%%%%%%%WRITING THE DATA TO A TEXT FILE%%%%%%%%%%%%%%
if  0
%     fid = fopen('ParameterValues.txt','wt');
%     fprintf(fid,'%f \n',Matrix);
%     fclose(fid);
writematrix(Matrix,'ParameterValues.txt','Delimiter',' ')
writematrix(Matrix2,'ParameterValues2.txt','Delimiter',' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
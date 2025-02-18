a1 = 1e-3;

R_vals = 10.^[0.1, 0.2, 0.5, 1, 2,3,4,5,6,7,8];

Int_values = zeros(length(R_vals),1);
for i=1:length(R_vals)
    R = R_vals(i);
    dmin = @(q) -sqrt(R*R-q.*q);
    dmax = @(q) sqrt(R*R-q.*q);
    qmin = -R;
    qmax = R;
    fun4 = @(q,d) (d.*(q.^2 + d.^2)./ sqrt(4*a1^2 * (q.^2 + d.^2) + 1));
    int4 = integral2(fun4,qmin,qmax,dmin,dmax,'RelTol',1e-5);
    
    Int_values(i) = int4;
end
figure;
loglog(R_vals,Int_values);
grid on;
grid minor;
ylabel('Integral values');
xlabel('R values');
figure;
fsurf(fun4, [-1 1 -50 50])
%int4 = integral2(fun4,qmin,qmax,dmin,dmax,'RelTol',1e-1,'Method','Iterated')
% int4_1 = integral2(fun4,qmin,qmax,dmin,dmax,'RelTol',1e-5)
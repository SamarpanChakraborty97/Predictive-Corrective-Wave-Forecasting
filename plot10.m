function h = plot10(A,plottingParams)

% color scheme
% 0 is not accepted
% 1 means use xi for the color
% 2-6 black - white

%determine radus of particles



%columns of A
%[x y vx vy rho color]

A = A(A(:,1)>plottingParams.axis(1),:);
A = A(A(:,1)<plottingParams.axis(2),:);
A = A(A(:,2)>plottingParams.axis(3),:);
A = A(A(:,2)<plottingParams.axis(4),:);

nPts  = 32;
N     = length(A(:,1));
h     = zeros(N,1);
theta = linspace(0,2*pi,nPts)';

%assemble polygons into a singel vector and plot it all at once
x = zeros(nPts,N);
y = zeros(nPts,N);

for ind1 = 1:N
    px = A(ind1,1);
    py = A(ind1,2);
    r = plottingParams.h;
    
   x(:,ind1) =  r*cos(theta)+px;
   y(:,ind1) =  r*sin(theta)+py;
   xi = A(:,5)'; %1D stress -> force/radius
    
end


C = A(:,6);

rgbColor = [0 .25 .5 .75 1];

%plot xi
fill(x(:,C==7),y(:,C==7),'b','edgecolor','none');
hold on


%plot gray1
fill(x(:,C==2),y(:,C==2),rgbColor([1 1 1]));

fill(x(:,C==3),y(:,C==3),rgbColor([2 2 2]));
fill(x(:,C==4),y(:,C==4),rgbColor([3 3 3]));
fill(x(:,C==5),y(:,C==5),rgbColor([4 4 4]));
fill(x(:,C==6),y(:,C==6),rgbColor([5 5 5]));

%plot gray2

    
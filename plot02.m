function h = plot02(A,plottingAxis,varargin)

%color scheme
% 0 is not accepted
% 1 means use xi for the color
% 2-6 black - white

%columns of A
%[x y vx vy rho color]

%filter down A
%A = A(A(:,1)>plottingAxis(1),:);
%A = A(A(:,1)<plottingAxis(2),:);
%A = A(A(:,2)>plottingAxis(3),:);
%A = A(A(:,2)<plottingAxis(4),:);

plottingSymbol = '.';
if length(varargin)==1
    plottingSymbol = varargin{1};
end

N     = length(A(:,1));
C     = A(:,4);

rgbColor = [0 .25 .5 .75 1];

x = A(:,1);
y = A(:,2);

%plot the free particles
plot(x(C==7),y(C==7),['b' plottingSymbol]);
plot(x(C==7), y(C==7), '.b','MarkerSize',3);
hold on


%plot gray1
plot(x(C==2),y(C==2),plottingSymbol,'color',rgbColor([1 1 1]));

plot(x(C==3),y(C==3),plottingSymbol,'color',rgbColor([2 2 2]));
plot(x(C==4),y(C==4),plottingSymbol,'color',rgbColor([3 3 3]));
plot(x(C==5),y(C==5),plottingSymbol,'color',rgbColor([4 4 4]));
plot(x(C==6),y(C==6),plottingSymbol,'color',rgbColor([5 5 5]));
%plot gray2
h = 0; %for backwards compatability
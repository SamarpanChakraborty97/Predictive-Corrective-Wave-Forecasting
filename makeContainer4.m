%make container1 makes a conical container with an opening a the bottom



function s = makeContainer4(s);

%tempWallx = s.free.xSpan(2);

%make a fixed container
% H = 0.1;
% H_by_dx = 20;
% dx = H/H_by_dx;
dx = 0.01;
%h = 2 * dx;

%length1 = 25;
length1 = 50;
%length2 = 4;
%depth = 0.4;

dxc = 0.3*dx;
h = 1.5 * dx;
yval = 1.0;

%y-value of wall 1
y1 = [yval:-dxc:0]';
%x-values of wall 1
x1 = -0.8 * ones(length(y1),1);

x2 = [-0.8+dxc:dxc:length1+h]';
%x2 = [dx:dx:0.5-dx]';
y2 = zeros(length(x2),1);

y3 = [yval:-dxc:0]';
%x-values of wall 1
x3 = (length1+h) * ones(length(y3),1);

y4 = [yval:-dxc:0]';
x4 = zeros(length(y4),1);

% 
% y5 = [1.2-dxc:-dxc:dxc]';
% x5 = (tempWallx+2*dx)*ones(length(y5),1);
% %now put them all together


pos = [[x1 y1];[x2 y2];[x3 y3];[x4 y4]];
pos2 = [[x1 y1];[x2 y2];[x3 y3]];

% r2=find(pos(:,1)==tempWallx+2*dx & pos(:,2)==dxc)
% r1 = find(pos(:,1)==tempWallx+2*dx & pos(:,2)==1.2-dxc)



nC = length(pos(:,1));

Piston_particle_start = length(pos2(:,1))+1
Piston_particle_end = length(pos(:,1))

s.comp.nConstrainedParticles = nC;

s.containerParticles.pos   = pos;
s.containerParticles.vel   = zeros(nC,3);
%s.containerParticles.mass  = s.free.mass*ones(nC,1);   %mass                         [kg/m^2]
%s.containerParticles.smoothingLength = s.free.smoothingLength*ones(nC,1);  %smoothingLength - mean, variance [m]
s.containerParticles.color = s.cons.color*ones(nC,1);


% plot(pos(:,1),pos(:,2),'k.')
% hold on 
% axis equal
% axis([-0.2 3.5 -0.1 2.1])

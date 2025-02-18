%make container1 makes a conical container with an opening a the bottom
function s = makeContainer3_load(s);

load ('PosContPart4.mat')
pos = Contparticles;


nC = length(pos(:,1));

s.comp.nConstrainedParticles = nC;

s.containerParticles.pos   = pos;
s.containerParticles.vel   = zeros(nC,3);
%s.containerParticles.mass  = s.free.mass*ones(nC,1);   %mass                         [kg/m^2]
%s.containerParticles.smoothingLength = s.free.smoothingLength*ones(nC,1);  %smoothingLength - mean, variance [m]
s.containerParticles.color = s.cons.color*ones(nC,1);


% plot(pos(:,1),pos(:,2),'k.')
% hold on 
% axis equal
% axis([-2 6 -2 2])
 









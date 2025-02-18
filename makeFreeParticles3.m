%function to assign the locations of free particles

function s = makeFreeParticles3(s)



%compute free particle positions
%nP = s.free.nSpan(1)*s.free.nSpan(2);

% H = 0.66;
% H_by_dx = 20;
% dx = H/H_by_dx;
% h = 1.5 * dx;

% length1 = 100;
% depth = 23;

x = linspace(s.free.xSpan1(1),s.free.xSpan1(2),s.free.nSpan1(1));
y = linspace(s.free.ySpan1(1),s.free.ySpan1(2),s.free.nSpan1(2));

[X,Y] = meshgrid(x,y);
U = reshape(X,[],1);
V = reshape(Y,[],1);
pos1 = [U V];
%size(pos1)

% % %%%TRIANGULAR PART OF THE NUMERICAL TANK%%%
% % length2 = 4;
% % x2 = linspace(s.free.xSpan2(1),s.free.xSpan2(2),s.free.nSpan2(1));
% % y2 = linspace(s.free.ySpan2(1),s.free.ySpan2(2),s.free.nSpan2(2));
% % 
% % [X2,Y2] = meshgrid(x2,y2);
% % 
% % A = reshape(X2,[],1);
% % B = reshape(Y2,[],1);
% % 
% % pos2 = [A B];
% % x1 = length1;
% % y1 = h;
% % x2 = length1+ length2;
% % y2 = depth;
% % m = (y2 - y1)/(x2 - x1);
% % 
% % pos2((pos2(:,2)- y1) < m*(pos2(:,1)- x1),:)=[];
% % %size(pos2)
% % %%%TRIANGULAR PART MAKING ENDS HERE%%%

position = [pos1];
s.freeParticles.pos = position;
s.freeParticles.nP = length(position);
nP = s.freeParticles.nP;
s.freeParticles.color = s.free.color*ones(length(X(:)),1);

%compute free particle characteristics
if length(s.free.smoothingLength)==1
    s.freeParticles.smoothingLength = ones(nP,1)*s.free.smoothingLength;
elseif length(s.free.smoothingLength)==2
    standardDev        = sqrt(s.free.smoothingLength(2));
    %s.freeParticles.ro = s.free.ro(1)+standardDev*randn(nP,1);
    %s.freeParticles.ro = 
    
    ro = s.free.smoothingLength(1)+standardDev*randn(nP,1); %Np x 1
    %if any particles are smaller than 2 std dev made them 2 std dev
    lowerCutoff = s.free.smoothingLength(1)-2*standardDev;
    upperCutoff = s.free.smoothingLength(1)+2*standardDev;
    ro(ro<lowerCutoff) = lowerCutoff;
    ro(ro>upperCutoff) = upperCutoff;
    s.freeParticles.smoothingLength = smoothingLength;

end

s.freeParticles.mass    = s.free.mass;
s.freeParticles.smoothingLength = s.free.smoothingLength;
s.freeParticles.cSound  = s.free.cSound;
s.freeParticles.nu      = s.free.nu;    
s.freeParticles.delta   = s.free.delta;    
s.freeParticles.muAlpha = s.free.muAlpha;     
s.freeParticles.muBeta  = s.free.muBeta;     %viscosity beta
s.freeParticles.epsilon = s.free.epsilon;     %epsilon XSPH
s.freeParticles.rRef    = s.free.rRef; %reference
%s.freeParticles.mass    = s.free.mass*ones(s.freeParticles.nP,1); %mass

%keyboard
% plot(s.freeParticles.pos(:,1),s.freeParticles.pos(:,2),'k.')
% axis equal
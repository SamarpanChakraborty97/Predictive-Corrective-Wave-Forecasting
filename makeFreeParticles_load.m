%function to assign the locations of free particles

function s = makeFreeParticles_load(s)


load ('PosFreePart4.mat')

X       = Freeparticles(:,1);
Y       = Freeparticles(:,2);
%compute free particle positions
nP = length(X);

% x = linspace(s.free.xSpan(1),s.free.xSpan(2),s.free.nSpan(1));
% y = linspace(s.free.ySpan(1),s.free.ySpan(end),s.free.nSpan(2));
% 
% [X,Y] = meshgrid(x,y);


s.freeParticles.pos = [X(:) Y(:)];
s.freeParticles.nP = length(X(:));

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

s.freeParticles.cSound  = s.free.cSound;     
s.freeParticles.nu      = s.free.nu;    
s.freeParticles.delta   = s.free.delta;    
s.freeParticles.muAlpha = s.free.muAlpha;     
s.freeParticles.muBeta  = s.free.muBeta;     %viscosity beta
s.freeParticles.epsilon = s.free.epsilon;     %epsilon XSPH
s.freeParticles.rRef    = s.free.rRef; %reference
s.freeParticles.mass    = s.free.mass;
s.freeParticles.smoothingLength = s.free.smoothingLength;

%keyboard
%plot(s.freeParticles.pos(:,1),s.freeParticles.pos(:,2),'k.')
%axis equal




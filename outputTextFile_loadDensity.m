function outputTextFile_loadDensity(s)

%output # free particles
%free particle position X
%free particle position Y
%free particle mass
%free particle radius
%free parameters

%output # constrained particles
%constrained particle position X
%constrained particle position Y
%constrained particle
%constrained particle radius
%constrained parameters
%constrained particles movement

%computational parameters




cd(s.path.pathToCase);
FID = fopen(s.path.nameOfInputFile,'w');
%FID = 1;

%free particles
nFree = length(s.freeParticles.pos(:,1));
fprintf(FID,'%d\n',-2);                               %tag to load density
fprintf(FID,'%10.0u\n',nFree);                        %#free
fprintf(FID,'%16.16f\n',s.freeParticles.pos(:,1));    %free x-pos
fprintf(FID,'%16.16f\n',s.freeParticles.pos(:,2));    %free y-pos
fprintf(FID,'%16.16f\n',s.freeParticles.density(:));     %free denisty
%fprintf(FID,'%16.16f\n',s.freeParticles.mass(:));     %free radius
%fprintf(FID,'%16.16f\n',s.freeParticles.smoothingLength(:));    %free radius
fprintf(FID,'%10.0u\n',7*ones(nFree,1));                      %free color

fprintf(FID,'%16.16f\n',s.freeParticles.mass);        % particle mass
fprintf(FID,'%16.16f\n',s.freeParticles.smoothingLength);        % particle smoothing length
fprintf(FID,'%16.16f\n',s.freeParticles.cSound);          %numerical speed of sound
fprintf(FID,'%16.16f\n',s.freeParticles.nu);          %free particles nu
fprintf(FID,'%16.16f\n',s.freeParticles.delta);          %free particles delta-SPH term 
fprintf(FID,'%16.16f\n',s.freeParticles.muAlpha);          %free particles alpha viscosity
fprintf(FID,'%16.16f\n',s.freeParticles.muBeta);      %free particles beta viscosity
fprintf(FID,'%16.16f\n',s.freeParticles.epsilon);     %free particles epsilon
fprintf(FID,'%16.16f\n',s.freeParticles.rRef);        %reference density

%constrained particles
nContainer = length(s.containerParticles.pos(:,1));
fprintf(FID,'%10.0u\n',nContainer);                        %#container
fprintf(FID,'%16.16f\n',s.containerParticles.pos(:,1)); %container x-pos
fprintf(FID,'%16.16f\n',s.containerParticles.pos(:,2)); %container y-pos
fprintf(FID,'%16.16f\n',s.containerParticles.density(:));     %container denisty
%fprintf(FID,'%16.16f\n',s.containerParticles.mass(:));    %container mass
%fprintf(FID,'%16.16f\n',s.containerParticles.smoothingLength(:));    %container smoothingLengh
fprintf(FID,'%10.0u\n',2*ones(nContainer,1));    %constrained color
%fprintf(FID,'%10.0u\n',nFree+(1:nContainer)');    %constrained color

%locations to be measured
% nMeasure = length(s.measuredLocations.pos(:,1));
% fprintf(FID,'%10.0u\n',nMeasure);                        %measured points
% fprintf(FID,'%16.16f\n',s.measuredLocations.pos(:,1)); %measured locations x-pos
% fprintf(FID,'%16.16f\n',s.measuredLocations.pos(:,2)); %measured locations y-pos
% fprintf(FID,'%16.16f\n',s.measuredParticles.density(:));     %measured denisty
% fprintf(FID,'%16.16f\n',s.measuredLocations.mass(:));    %measured locations mass
% fprintf(FID,'%16.16f\n',s.measuredLocations.smoothingLength(:));    %measured locations smoothingLengh
% fprintf(FID,'%10.0u\n',5*ones(nMeasure,1));    %constrained color


%fprintf(FID,'%10u\n',3);                       %Xfunction 3 - p1*sin(2*pi*t*p2+p3)
% fprintf(FID,'%10u\n',10);
% %fprintf(FID,'%10u\n',1);                      %p1
% %fprintf(FID,'%10u\n',0.8);                      %p2
% %fprintf(FID,'%10u\n',0);                      %p3
% fprintf(FID,'%10u\n',1);                        %Yfunction 10 - do nothing
% fprintf(FID,'%10u\n',3);
% fprintf(FID,'%10u\n',0);
% %fprintf(FID,'%10u\n',2);
% %fprintf(FID,'%10u\n',312);
% fprintf(FID,'%10u\n',3683);                        %range1
% fprintf(FID,'%10u\n',4181);                      %range2

fprintf(FID,'%10u\n',3);                       %Xfunction 3 - p1*sin(2*pi*t*p2+p3)
fprintf(FID,'%16.16f\n',-5.0);            %wavemaker location
%fprintf(FID,'%16.16f\n',0.5);                  %frequency
%fprintf(FID,'%16.16f\n',0);                    %phase
%fprintf(FID,'%16.16f\n',4.52);                  %wavelength
%fprintf(FID,'%16.16f\n',0.66);                  %depth
%fprintf(FID,'%16.16f\n',1.02); 
fprintf(FID,'%10u\n',10);                      %Yfunction 10 - do nothing
fprintf(FID,'%10u\n',12880);                       %range1
fprintf(FID,'%10u\n',13890);                      %range2

%fprintf(FID,'%10u\n',10);                      %Xfunction 10 - do nothing - takes no parameters
%fprintf(FID,'%10u\n',10);                      %Yfunction 10 - do nothing
%fprintf(FID,'%10u\n',1);                       %range1
%fprintf(FID,'%10u\n',10);                      %range2
%

%computational parameters
fprintf(FID,'%10u\n',0);                      %no more functions
fprintf(FID,'%16.16f\n',s.comp.gravity);      %gravity
fprintf(FID,'%10.0u\n',s.comp.nT);            %# time steps
fprintf(FID,'%16.16f\n',s.comp.dt);           %dt
fprintf(FID,'%10.0u\n',s.comp.storageStride); %storage Stride

P = readmatrix('InitP.txt');
L = length(P(:,1));
%fprintf(FID,'%10.0u\n',L);                    %# number of frequency steps
for i=1:4
    for j = 1:L
        fprintf(FID,'%16.16f\n',P(j,i));
    end
end

fclose(FID);




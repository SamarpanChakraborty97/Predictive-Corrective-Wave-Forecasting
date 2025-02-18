close all
clear all

%%%%%%%%%%%%%%%%
%RELATIVE PATHS%
%%%%%%%%%%%%%%%%
pathToSPHCode = 'D:/CCSPH_GPU_GridStudy/';
%pathToSPHCode = '/~/scratch/SPH2D/';
%pathToSPHCode = '../';
pathToHere = mfilename('fullpath');
pathToDataFileStorage = 'dataFiles/';
nameOfInputFile       = 'SPHInputFile_Initialization2.txt';

%%%%%%%%%%%%%%%%
%FREE PARTICLES%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%PART1 OF NUMERICAL WAVE TANK%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H = 0.1;
%H_by_dx = 15;
%dx = H/H_by_dx;
H = 1.98;

H_dx = 10;
dx = H / H_dx;

%depth = 0.25;
depth = 27.4;

h = 1.0 * dx;
length1 = 300;
ch = sqrt(2 * 9.81 * depth);
vmax = 1.0;

s.free.xSpan1  = [h length1];
s.free.ySpan1  = [h depth];

xnSpan1 = length(h:dx:length1);
ynSpan1 = length(h:dx:depth);
s.free.nSpan1  = [xnSpan1 ynSpan1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%PART2 OF NUMERICAL WAVE TANK%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length2 = 100;
s.free.xSpan2 = [length1+dx length1+length2];
s.free.ySpan2 = [h depth];

xnSpan2 = length(length1+dx:dx:length1+length2);
ynSpan2 = length(dx:dx:depth);

s.free.nSpan2 = [xnSpan2 ynSpan2];

%define free particle positions
% s.free.xSpan  = [h width];
% s.free.ySpan  = [h height];
% xnSpan        = length(h:dx:width)
% ynSpan        = length(h:dx:height)
% %pause
% % xnSpan = 100;
% % ynSpan = 76; 
% s.free.nSpan  = [xnSpan ynSpan];

% nfreeParts= xnSpan*ynSpan;
% desiredHeight = 2*height;
% containerWidth = height;
% 
% desiredHeight = 0.025;
% containerWidth = 0.21
rho_a = 1000;
Va = dx ^ 2;
m = rho_a * Va;

%massParticle = (desiredHeight)*(containerWidth)*1000/nfreeParts
% 
% %parameters from book:
s.free.mass            = m;       %mass  [kg/m^2]                
s.free.smoothingLength = h;         
%s.free.vf          = 1.5;   %maximum fluid velocity
s.free.cSound      = 10 * max(ch,vmax);
s.free.nu          = 1 * 10^(-6);
s.free.delta       = 0.1; %kinematic viscosity
s.free.muAlpha     = 0; %-0.4
s.free.muBeta      = 1.2;  %muBeta -0.5
s.free.epsilon     = 0;  %epsilon in XSPH %0.25
s.free.rRef        = 1000; %reference density
s.free.color       = 1;                                   

%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINED PARTICLES%
%%%%%%%%%%%%%%%%%%%%%%%

%define container parameters
%container parameters are fixed
s.cons.color = 2;
s.measured.color = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTATIONAL PARAMETERS%
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.comp.gravity        = -9.8;

%export parameters
s.comp.storageStride  = 20000;

%time parameters
time_simulation       = 10.0;
s.comp.dt             = 5*10^(-5);
s.comp.nT             = ceil(time_simulation/s.comp.dt);


%%%%%%%%%%%%%%%%
%END USER INPUT%
%%%%%%%%%%%%%%%%

%create paths
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
s.path.pathToCase = ('.\');
s.path.nameOfInputFile  = nameOfInputFile;
pathToCCode       = [pathToSPHCode];

%make the data structures
s = makeFreeParticles2(s);
%s = makeFreeParticles_load(s);
%s = pressurePoint(s);
s = makeContainer3 (s);
%s = makeContainer3_load(s);
%s = LoadDensity(s); % problem with loading density when reverting back to makeContainer2 
plottingAxis = [-2 20 -0.1 2.1];
disp(s.freeParticles.nP)
length(s.containerParticles.pos(:,1))
% disp([s.freeParticles.nP+ length(s.containerParticles.pos(:,1))])

%2-312
%566- 695

if 1 and(1,length(s.containerParticles.pos(:,1))>853)
    ix = [2:312];
    plot(s.freeParticles.pos(:,1),s.freeParticles.pos(:,2),'b.')
    hold on
    plot(s.containerParticles.pos(:,1),s.containerParticles.pos(:,2),'k.')
    %plot(s.measuredLocations.pos(:,1), s.measuredLocations.pos(:,2),'r.')
    %plot(s.containerParticles.pos(ix,1),s.containerParticles.pos(ix,2),'ro')
    axis equal
end

% return

%ouput text file
outputTextFile(s)
%outputTextFile_loadDensity(s)
%outputTextFile_loadDensity2(s)

%execute C code
I = regexp(pathToCCode,'/'); %replace / with \
pathToCCode(I) = '\'; 

I = regexp(pathToDataFileStorage,'/'); %replace / with \
pathToDataFileStorage(I) = '\'; 

%dos(['cd ' pathToCCode '&&Trial_WCSPH_cpu.exe' pathToHere nameOfInputFile ' ' pathToHere pathToDataFileStorage ' &'],'-echo');
%dos(['cd ' pathToCCode '&&Trial_WCSPH_cpu.exe ../inputFile.txt']);
%dos_command = ['"C:/Users/schakr18/source/repos/Trial_WCSPH_cpu/Debug/Trial_WCSPH_cpu.exe" "C:/Users/schakr18/source/repos/Trial_WCSPH_cpu/Debug/SPHInputFile_Initialization.txt"'];
%system(dos_command);
%system(['cd ' pathToCCode '&&Trial_WCSPH_cpu.exe' pathToHere nameOfInputFile]);

function s = LoadDensity(s)
load ('DensityFreePart4.mat')
s.freeParticles.density = DensityFreeParticles(:);

load ('DensityContPart4.mat')
s.containerParticles.density = DensityContainerParticles(:);
% 
% load('DensityMeasPart4.mat')
% s.measuredParticles.density = DensityMeasuredParticles(:);
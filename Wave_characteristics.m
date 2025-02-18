clc
clear all;

%%%%%%%%% WAVE CHARACTERISTICS %%%%%%%%%%
unit_conv = 3.28084;
T = 2;
g = 9.81 * unit_conv;
gT2= g/(T^2);

% req_varH = 0.00012;
% req_varD = 0.0012;
% 
% H = gT2 * req_varH;
% d = gT2 * req_varD;


H = 0.1 * unit_conv;
H/(T^2)

d = 0.5 * unit_conv;
d/(T^2)

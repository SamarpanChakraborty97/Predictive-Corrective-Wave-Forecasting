clear all
 h = 0.1;
 L = 0.6;
 k = 2*pi/L;
 kh = k*h
 S = 0.04;
 ratio = 2*(cosh(2*kh)-1)/(sinh(2*kh)+2*kh);
 H = ratio*S
 
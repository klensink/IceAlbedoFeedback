function [OLR] = OLR(T,T0,CO2,CO20,alpha)
 
%T0 is ref T in absolute T
%C_C0 is pCO2 change relative to baseline
%S0 = solar constant
%alpha = albedo
 S0=1368; 
A=(S0/4)*(1-alpha)-4.*log2(CO2./CO20);
OLR=A+2.*(T-T0);

global c e me B theta lnD eps0 mu0 re Z;
global r_torus R_torus rcameraExt rcameraInt;
global Lp0 Mpt0 Mpv0;
global windowSize cFiltro dFiltro;
global scaledFactorElectronNumber scaledFactorVloop scaledFactorNeutrnfc144P scaledFactorNeutrnfc144M;
%% 
c=3*10^8;               %speed of light
e=-1.602*10^(-19);      %electron charge
me=9.109*10^(-31);       %
eps0=8.854*10^(-12);    %
mu0 = 4*pi*10^(-7);     %Permeabilità magnetica H/m,
re= 2.8179*10^(-15); %e^2/(4*pi*eps0^2*me*c^2);
Z=1;

%%TOKAMAK
B=5.5;                  %toroidal magnetic field
theta=0.5;              %pitch angle 30° 
lnD=15;                 %Coulomb logarithm   ??
Lp0=2*10^(-6);
Mpt0=110*10^(-6);
Mpv0=80*10^(-6);
R_torus=0.935;
r_torus=0.305;
rcameraExt=0.935+0.305;
rcameraInt=0.935-0.305;

%%SIMULATION
windowSize=6;
cFiltro = 45; %try also c=1
dFiltro = 90; %d = 2 and see the difference

%% 
scaledFactorElectronNumber=1e-14;
scaledFactorVloop=1e4;
scaledFactorNeutrnfc144M=0.54e-9;
scaledFactorNeutrnfc144P=0.54e5;

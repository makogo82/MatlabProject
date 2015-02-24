global c e me B theta lnD eps0 mu0 re Z;
global Lp0 Mpt0 Mpv0;
%% 
c=2.99792*10^8;               %speed of light
e=-1.602*10^(-19);      %electron charge
me=9.109*10^(-31);       %
eps0=8.854*10^(-12);    %
mu0 = 4*pi*10^(-7);     %Permeabilità magnetica H/m,
re= 2.8179*10^(-15); %e^2/(4*pi*eps0^2*me*c^2);
Z=1;

%%TOKAMAK
B=5.5;                  %toroidal magnetic field
theta=deg2rad(30)       %pitch angle 30° 
lnD=10;                 %Coulomb logarithm   ??
Lp0=2*10^(-6);
Mpt0=110*10^(-6);
Mpv0=80*10^(-6);


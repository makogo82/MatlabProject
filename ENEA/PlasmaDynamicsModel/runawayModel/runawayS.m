function [dx,y] = runawayS(t,x,u,p1,p2,varargin)
%% INPUT 
% simple model EllN TeMaxN neMaxN
ne=u(1);
Ell=u(2);
%dynamic model+
dx = x*(-p1*ne  + p2*Ell);
%% Output [position, current]
y = x;


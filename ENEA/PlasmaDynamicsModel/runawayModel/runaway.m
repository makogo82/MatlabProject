function [dx,y] = runaway(t,x,u,p1,p2,p3,varargin)
%% INPUT 
% simple model EllN TeMaxN neMaxN
ne=u(1);
Ell=u(2);
%dynamic model+
dx = -p1*ne  + p2*Ell + p3*x;
%% Output [position, current]
y = x;


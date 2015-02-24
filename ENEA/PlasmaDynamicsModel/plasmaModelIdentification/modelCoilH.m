function [dx,y] = modelCoilH(t,x,u,p1,p2,p3,p4,p5,p6,varargin)
%% INPUT 

IT=u(1);
IV=u(2);
IF=u(3);
IH=u(4);
%dynamic model
dx1 = -p1*x(1) + p4*IH + x(2);
dx2 = -p2*x(1) + p5*IH + x(3);
dx3 = -p3*x(1) + p6*IH; 


dx=[dx1;dx2;dx3];
%% Output [position, current]
y = x(1);


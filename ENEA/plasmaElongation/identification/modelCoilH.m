function [dx,y] = modelCoilH(t,x,u,p1,p2,p3,p4,varargin)

%% INPUT

zs2 = u(5);
zs1 = u(4);
IF = u(3);
IH = u(2);
IP = u(1);

% dynamic model 
% tengo conto dell'effetto fi F
dx1 = x(2);

b_pl = (zs1-zs2)/2;

xV = x(1) + b_pl;

xV = x(1);

%EFFETTO DOVUTO A IH
f_IH_IP = p1*IP*IH*( 1 / (p2-xV) - 1 / (-p2-xV) );

%EFFETTO DOVUTO A IF
f_IF_IP = p3*IP*IF*( - 1 / (p4+x(1)) + 1 / (p4-x(1)) ); 

dx2 = f_IH_IP + f_IF_IP;

dx=[dx1;dx2];
%% Output [position, current]
y = x(1);


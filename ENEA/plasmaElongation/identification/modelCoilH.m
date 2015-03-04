function [dx,y] = modelCoilH(t,x,u,p1,p2,p3,p4,varargin)

%% INPUT

zs1 = u(4);
zs2 = u(5);
dez = u(3);
IH = u(2);
IP = u(1);

% dynamic model 
% tengo conto dell'effetto fi F
dx1 = x(2);
dx2 = p1*IH*IP*(1/(0.8+x(1))+1/(0.8-x(1)))+ p2*(IV/(0.8-x(1)))^2;

p1 = 1E-5;
p2 = c_2*0.001; 
p3 = 0.4;
p4 = 0.8;
b_pl = (zs1-zs2)/2;

f_IH_IP = p1*IP*IH*( 1 / (p2-(dez+b_pl)) - 1 / (-p2-(dez+b_pl)) );
f_IF_IP = -p3*IP*IF*( 1 /(p4-(dez+b_pl)) + 1 / (-p4-(dez+b_pl)) )';

dx2 = f_IH_IP + f_IF_IP;

dx=[dx1;dx2];
%% Output [position, current]
y = x(1);


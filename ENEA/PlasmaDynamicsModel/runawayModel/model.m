function dy=model( t,y, newTime,ne, Ell)
dy = zeros(1,1);
p0=130;
p1=1.1*0.275/0.1843;

neI=interp1(newTime,ne,t);
EllI=interp1(newTime,Ell,t);

%[t neI EllI]

dy(1) = p0*y(1)*(-p1*neI*2E-21+EllI)


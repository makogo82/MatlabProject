%% About 37714 data
%- time is given in seconds, 
%- 'neMax' is the central electron density in cubic metres, 
%- 'n' is just a plotting aid (neMax*10^(-19)), 
%- 'TeMax' (central electron temperature) is in keV and 
%- variables 'a', 'b', 'c' and 'd' are auxiliary, used to 
% calculate 'Spri' (as a*b*c*d). 
% The values of 'neMax' are experimental from 0.3s, the rest 
% is extrapolation by the program. Same goes for 'Zeff' after 1s and 'TeMax' after 1.45s.

data.a = a;
data.b = b;
data.c = c;
data.d = d;
data.Ed = ED;
data.Ell= EllED;
data.n=n;
data.neMax=neMax;
data.Spri=Spri;
data.TeMax=TeMax;
data.time=time;
data.Zeff=Zeff;



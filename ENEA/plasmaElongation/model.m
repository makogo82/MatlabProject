%modello approssimato della posizione verticale
clear all; clc; close all;
c_2 = 1E-5;
c_2H = 1.7E-5;
c_2L = 0.3E-5;
c_2 = 1E-5;
c_3 = 0.8;
Ts = 0.0005;

shots =         [39091   39092  39094  38629  38628    38627  38404  38616 ];
shots_time =    [0.8    1.08    0.85   0.9025   0.5     1.14    0.8    0.7 ];
shots_endtime = [0.904   1.14   1      0.985   1.195    1.22   1.178   0.9 ];

shot =38629;
Ts = 0.5E-3;
samplingTime = Ts;
t_start = shots_time(find(shots==shot));
t_end = shots_endtime(find(shots==shot));
data = CleanFTUdata_v6(t_start,t_end,Ts,shot,1);
t=data.time;

c_dez = 2;
d_dez = 9;
y = data.dez/((1.5e-4) * (-1.0));
dy = mypseudo_derivative(t,y,c_dez,d_dez);



x1(1)=y(1);
x1L(1) =y(1);
x1H(1) =y(1);
x2(1) = dy(1);
x2L(1) =dy(1);
x2H(1) =dy(1);
%%
time(1) = t(1);
for h=1:length(t)-1
    IPLmis = data.IPLmis(h);
    Hmis = data.Hmis(h);
    Fmis = data.Fmis(h);
    t1 = c_2*IPLmis*Hmis*(1/(c_3-1.5*(x1(h))) - 1/(-c_3-1.5*(x1(h))) );
    t1L = c_2L*IPLmis*Hmis*(1/(c_3-1.5*(x1L(h))) - 1/(-c_3-1.5*(x1L(h))) );
    t1H = c_2H*IPLmis*Hmis*(1/(c_3-1.5*(x1H(h))) - 1/(-c_3-1.5*(x1H(h))) );
    c_4 = c_2*0.001;
    c_5 = 0.4; 
    t2 = -c_4*Fmis*IPLmis*(2*( 1/(c_5-1.5*(x1(h))) + 1/(-c_5-1.5*(x1(h)))))^2';
    t2L = -c_4*Fmis*IPLmis*(2*( 1/(c_5-1.5*(x1L(h))) + 1/(-c_5-1.5*(x1L(h)))))^2';
    t2H = -c_4*Fmis*IPLmis*(2*( 1/(c_5-1.5*(x1H(h))) + 1/(-c_5-1.5*(x1H(h)))))^2';
    %%

    % EULERO (funzioni)
    x2(h+1) = x2(h) + Ts*(t1+t2);
    x2L(h+1) = x2L(h) + Ts*(t1L+t2L);
    x2H(h+1) = x2H(h) + Ts*(t1H+t2H);
    x1(h+1) = x1(h) + Ts*x2(h);
    x1L(h+1) = x1L(h) + Ts*x2L(h);
    x1H(h+1) = x1H(h) + Ts*x2H(h);
    time(h+1) = time(h)+Ts;

end
%%
close all
figure('Name','velocity/position')
subplot(2,1,1)
plot(time,x2,time,x2L,'k-.',time,x2H,'r-.');
legend('x2','x2L','x2H','ddez')
grid on;
subplot(2,1,2)
plot(time,x1,time,x1L,'k-.',time,x1H,'r-.');
grid on;
legend('x1','x1L','x1H','dez')
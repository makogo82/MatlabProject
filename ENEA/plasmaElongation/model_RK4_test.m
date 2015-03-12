%modello approssimato della posizione verticale
clear all; clc; close all;
c_2H = 1.5E-5;
c_2L = 0.5E-5;
c_2 =  1E-5;
c_3 = 0.7;
c_4 = c_2*0.001;
c_5 = 0.4
Ts = 0.0005;

shots =         [39091   39092   39094  38629      38628    38627  38404  38616 ];
shots_time =    [0.8     1.12    0.9    0.82        0.85     1.15    0.8    0.7 ];
shots_endtime = [0.904   1.14    0.95   0.984      1.195       1.22   1.178   0.9 ];

shot =38627;

samplingTime = Ts;
t_start = shots_time(find(shots==shot));
t_end = shots_endtime(find(shots==shot));

dataAll = CleanFTUdata_v6(0.2,1.2,Ts,shot,1);
t = [t_start:Ts:t_end];

data.IPLmis = interp1(dataAll.time,dataAll.IPLmis, t);
data.Hmis   = interp1(dataAll.time,dataAll.Hmis,   t);
data.Fmis   = interp1(dataAll.time,dataAll.Fmis,   t);
data.dez    = interp1(dataAll.time,dataAll.dez,    t);
data.d_dez  = interp1(dataAll.time,dataAll.d_dez,    t);
data.dd_dez = interp1(dataAll.time,dataAll.dd_dez,    t);

x1o(1)=0;
x2o(1)  =  0;
x1E(1)=0;
x2E(1)  =  0;
%%
time(1) = 0;
Ts = 0.5E-2;
for h=1:2e5
    
    
    
    [x1o(h+1), x2o(h+1)]  = verticalModel_RK(time(h),Ts , x1o(h) , x2o(h) , 1,1,1, c_2, c_3 , c_4 , c_5  );
    [x1E(h+1), x2E(h+1)]  = verticalModel_Eulero(time(h),Ts , x1E(h) , x2E(h) , 1,1,1, c_2, c_3 , c_4 , c_5  );

    %RUNGE-KUTTA

    time(h+1) = time(h)+Ts;

end

%%
close all
figure('Name','velocity/position')
subplot(2,1,1)
plot(time,x2o,'r',time,x2E,'b-.');
legend('RK4','Eulero')
title('$\dot{dez}$','Interpreter','latex')
grid on;
subplot(2,1,2)
plot(time,x1o,'r',time,x1E,'b-.');
grid on;






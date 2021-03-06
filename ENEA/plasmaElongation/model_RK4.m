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
Ts = 0.5E-3;
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

x1(1)=data.dez(1);
x1o(1)=data.dez(1);
x1L(1) = data.dez(1);
x1Lo(1) = data.dez(1);
x1H(1) = data.dez(1);
x1Ho(1) = data.dez(1);
x2o(1)  =  data.d_dez(1);
x2(1)  =  data.d_dez(1);
x2L(1) = data.d_dez(1);
x2Lo(1) = data.d_dez(1);
x2H(1) = data.d_dez(1);
x2Ho(1) = data.d_dez(1);
%%
time(1) = t(1);
for h=1:length(t)-1
    
    IPLmis = data.IPLmis(h);
    Hmis = data.Hmis(h);
    Fmis = data.Fmis(h);
    
    [x, t1, t2 ] =    verticalModel(time(h),x1(h),IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5); 
    [xL, t1L, t2L ] = verticalModel(time(h),x1L(h),IPLmis,Hmis,Fmis, c_2L , c_3 , c_4 , c_5);
    [xH, t1H, t2H ] = verticalModel(time(h),x1H(h),IPLmis,Hmis,Fmis, c_2H , c_3 , c_4 , c_5);
        
    
    
    [x1o(h+1), x2o(h+1)]  = verticalModel_RK(time(h),Ts , x1o(h) , x2o(h) , IPLmis,Hmis,Fmis, c_2, c_3 , c_4 , c_5  );
    [x1Lo(h+1), x2Lo(h+1)]  = verticalModel_RK(time(h),Ts , x1Lo(h) ,   x2Lo(h) , IPLmis,Hmis,Fmis, c_2L, c_3 , c_4 , c_5  );
    [x1Ho(h+1), x2Ho(h+1)]  = verticalModel_RK(time(h),Ts , x1Ho(h) ,   x2Ho(h) , IPLmis,Hmis,Fmis, c_2H, c_3 , c_4 , c_5  );

    
    %RUNGE-KUTTA
      
    
    % EULER
    %velocity
    x2(h+1) =  x2(h) +  Ts*(t1+t2);
    x2L(h+1) = x2L(h) + Ts*(t1L+t2L);
    x2H(h+1) = x2H(h) + Ts*(t1H+t2H);
    %position
    x1(h+1) = x1(h)   + Ts*x2(h);
    x1L(h+1) = x1L(h) + Ts*x2L(h);
    x1H(h+1) = x1H(h) + Ts*x2H(h);
    %
    time(h+1) = time(h)+Ts;

end

%%
close all
figure('Name','velocity/position')
ax(1) = subplot(2,1,1)
plot(time,x2,'g',   time,x2o,'g:',...
     time,x2L,'k',time,x2H,'r',...
     time,x2Lo,'k:',time,x2Ho,'r-.',...
     time,data.d_dez,'c-.');
ylim([-11 11]);
xlim([t_start t_end]);
legend('x2Euler','x2RK4','x2LEuler','x2HEuler','x2LRK4','x2HRK4','ddez')
title('$\dot{dez}$','Interpreter','latex')
grid on;
plotParam;
ax(2) = subplot(2,1,2)
plot(time,x1,'g',time,x1L,'k',time,x1H,'r',...
     time,x1o,'g:',time,x1Lo,'k:',time,x1Ho,'r:',...
     time,data.dez,'c-.');
ylim([-0.2 0.2]);
xlim([t_start t_end]);
grid on;
legend('x1Euler','x1LEuler','x1HEuler','x1RK4','x1LRK4','x1HRK4','dez')
title('$dez$','Interpreter','latex')
linkaxes(ax,'x');
plotParam;






clc
clear all;
close all;
globalVar
t_i=0;
t_f=1;
Ts=0.0005;
local=0;
shot = 36634;
windowFilter=10
b=ones(1,windowFilter)/windowFilter;
dato2 = load(strcat('shot_', num2str(shot),'.mat'));  %dati zana
dato = dato2.Data;
%Data = CleanFTUdata_v6(t_i, t_f, Ts, shot, local) 
%plot(dato.IplOldMis.y);
time = [0.01:Ts:1.5]
Ipl=  interp1(dato.IPLmis.x,dato.IPLmis.y,time);

Vloop =  interp1(dato.oldVloop.x,dato.oldVloop.y,time);
ALVmis = interp1(dato.ALVmis.x,dato.ALVmis.y,time);
ALTmis = interp1(dato.ALTmis.x,dato.ALTmis.y,time);

%plot(dato.ALVmis.y); 
figure('Name','IPL');
plot(time,Ipl);
plotParam;
figure('Name','Vloop')
plot(time,Vloop);
plotParam;
%%

Rp=5e-6;
Mpt0=107E-6;
Mpv0=75E-6;

%%

figure('Name',strcat('ALVmis_ALTmis', num2str(shot)))
ax1(1) = subplot(1,1,1)
IplFilter = filtfilt(b,1,Ipl);
VloopFilter = filtfilt(b,1,Vloop);
RpFilter = filtfilt(b,1,VloopFilter./IplFilter);
plot(time,Vloop./Ipl, time, VloopFilter./IplFilter); 
legend('Vloop/IP','VloopFilter/IPFilter')
plotParam;

%%
figure('Name',strcat('DerivataALVmis_ALTmis', num2str(shot)))
ax2(1) = subplot(2,1,1)
ALVmisFilter = filtfilt(b,1,ALVmis); 
plot(time, ALVmis, time, ALVmisFilter); 
legend('dAlv')
plotParam;
ax2(2) = subplot(2,1,2)

ALTmisFilter = filtfilt(b,1,ALTmis); 
plot(time, ALTmis , 'r-.', time, ALTmisFilter); 
legend('ALTmis','ALTmisFilter')
plotParam;
linkaxes(ax2,'x');
%%
simObject =sim('plasmaModel',time);
figure('Name',strcat('IP_dIT,d_IV.shot_', num2str(shot)))
ax(1) = subplot(2,1,1)
plot(time,RpFilter); 
legend('R_P')
plotParam;
ax(2) =subplot(2,1,2)
plot(dIVFilter.Time, dIVFilter.Data, dITFilter.Time, dITFilter.Data);
legend('dIV','dIT')
plotParam;
linkaxes(ax,'x');
%%
thetaT = 25.533; 
%sparo 36569
figure('Name',strcat('ThetaV', num2str(shot)))
plot(time,IplFilter,'r-.')
hold on;
plot(dITFilter.Time', -dITFilter.Data'*thetaT - 15*dIVFilter.Data','k')
plotParam;


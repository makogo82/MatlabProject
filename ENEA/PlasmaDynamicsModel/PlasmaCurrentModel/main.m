clc
clear all;
close all;
globalVar
t_i=0;
t_f=1;
Ts=0.0005;
local=0;
shot = 31855;
dato = load(strcat('shot_', num2str(shot),'.mat'));  %dati zana
%Data = CleanFTUdata_v6(t_i, t_f, Ts, shot, local) 
%plot(dato.IplOldMis.y);
time = [0.4:Ts:1.4]
Ipl=  interp1(dato.IplOldMis.x,dato.IplOldMis.y,time);

Vloop =  interp1(dato.Vloop.x,dato.Vloop.y,time);
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
figure('Name','Vloop/Ipl')
plot(time,Vloop./Ipl);
plotParam;
Rp=2.5e-6;
%%
figure('Name','ALVmis_ALTmis')
plot(time,ALTmis,time,ALVmis);

%%
[T,X,Y] =sim('plasmaModel',time);

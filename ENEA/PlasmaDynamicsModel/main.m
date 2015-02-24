clc;clear all;close all;
Ts = 0.0005; 

shot = 38570
t_i = 0;
t_f = 1.81;
Ts = 0.0005;
% get idddata
shotData = dataProcessing(t_i,t_f,Ts,shot, 3, 0.05);
%%
figure('Name',strcat(['IP_=',num2str(shot)]))
plot(shotData.time , shotData.IPLmis , 'r', shotData.time, shotData.IPLmisLowFilter, 'b');
legend('normal','low');
grid on; title(strcat(['I_{plasma}. Sparo=',num2str(shot)]));

figure('Name',strcat(['IH_=',num2str(shot)]))
plot(shotData.time , shotData.IHmis , 'r', shotData.time, shotData.IHmisLowFilter, 'b');
legend('normal','low');
grid on; title(strcat(['I_{H}. Sparo=',num2str(shot)]));

figure('Name',strcat(['IT_=',num2str(shot)]))
plot(shotData.time , shotData.ITmis , 'r', shotData.time, shotData.ITmisLowFilter, 'b');
legend('normal','low');
grid on; title(strcat(['I_{T}. Sparo=',num2str(shot)]));

figure('Name',strcat(['IF_=',num2str(shot)]))
plot(shotData.time , shotData.IFmis , 'r', shotData.time, shotData.IFmisLowFilter, 'b');
legend('normal','low');
grid on; title(strcat(['I_{F}. Sparo=',num2str(shot)]));



%% Magnetic 
figure('Name','VB01')
plot(shotData.time, shotData.vslow);

%%
z = shotData.iddDataLow; 
%identificazione del modello n4sid
sys = n4sid(z,6)
%figure('Name', 'Compare model whit data. Shot ');
compare(sys,z)


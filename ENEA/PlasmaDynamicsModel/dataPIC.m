close all;
clear all;
data = load('piclldata_command_T7_T8.mat')
datam=data.datam;
fine = 8e3;
t=datam(1:fine,1);
T1=datam(1:fine,2);
T2=datam(1:fine,3);
u=datam(1:fine,4);
figure(1)
ax(1)=subplot(2,1,1)
plot(t,u)
ax(2)=subplot(2,1,2)
plot(t,T1,t,T2)
linkaxes(ax)
%%
DAT = iddata(T2,u,t(2)-t(1))

%%
sys = n4sid(DAT,1);

compare(sys,DAT);

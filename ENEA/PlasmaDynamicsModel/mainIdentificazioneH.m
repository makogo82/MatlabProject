clc;clear all;close all;
Ts = 0.0005; 
shot = 38869
t_i = 0;
t_f = 1.81;
Ts = 0.0005;
% get idddata


tempus = load([ 'shot_' num2str(shot) '.mat']);
data = tempus.Data;

figure('Name',strcat(['HDEZ_=',num2str(shot)]))
subplot(2,1,1);
uh = data.ALHcalc.y;
plot(data.ALHcalc.x,uh);
legend('ALHcalc'); grid on;
subplot(2,1,2);
dez = data.DEZ.y
plot(data.DEZ.x,dez);
legend('DEZ'); grid on;

t=[0:0.0005:1.315];
newUh=interp1(data.ALHcalc.x,uh,t);
newDez=interp1(data.DEZ.x,dez,t);


figure('Name',strcat(['New HDEZ_=',num2str(shot)]))
subplot(2,1,1);
plot(t,newUh-newUh(1));
legend('ALHcalc'); grid on;
subplot(2,1,2);
plot(t,newDez);
legend('DEZ'); grid on;

%%
iddataH = iddata(newDez', newUh' , Ts, 'Name','identification');
plot(iddataH);

%%
%identificazione del modello n4sid
[sys,x0] = n4sid(iddataH,2)
%figure('Name', 'Compare model whit data. Shot ');
compare(sys,iddataH);

%%
close all;
[y,t,x] = lsim(sys,newUh-newUh(1),t); 
plot(t,iddataH.OutputData,'k', t,y, 'r')
legend('Measured', 'Simulated')

%%
% Kp = 1;
% Ki = 0;
% Kd = 1.91;
% C = pidstd(Kp,Ki,Kd,4,Ts)
% pidTF = tf(C);
% sysTF=feedback(pidTF*sys,1);
% impulse(sysTF)



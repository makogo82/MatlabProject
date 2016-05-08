

clc
clear all
close all
cols=['r' , 'm' ,'g' ,'y' ,'c' ,'b', 'k', 'r' , 'm'];

shot = [39903]%

% 38513
% 38519
% 38864
% 38865
% 38869
% 38870
% 39562
% 39903


shot = [37761];
j=1;
shot = [36767];
dataShot = ftudata(37761, '0_MACC.FBCOR-IP.A');
dataShot2 = ftudata(37761, 'MARTEFE.FBIP_M_I');
plot(dataShot.x, dataShot.y,dataShot2.x, dataShot2.y)

shot = [35965];
%shot = [38513];
shot = [37761];
%% no good for marte signals [33445  33485 34477 33659 33485 31630
delete signalsremote.mat;
for k = 1:length(shot)

    generate_datafile(shot(k));
    
    %%
    %cfg='../../MARTe-GenericTimerDrv-LcsmGAM-ExtSeek.cfg'
    %cfg ='../../MARTe-GenericTimerDrv-FTU-WithHModel.cfg'
    cfg = 'MARTe-Runaway-38513_PIDenanched38865.cfg';

    system(['./startSim_remote.sh ' cfg])
    
end
%%
%% show_figure;

cq_dIpl_threshold = -1E7;
TRIGGER_threshold_vloop = -4.5;
TRIGGER_threshold_ipl_vloop_activation =2E4;

load signalsremote.mat;
%%
close all;
figure(1)

time = CLK*1E-6;
t_i = 10.0;
t_f = 10.85;
signIPl = -1
dim = ones(1,length(time));
ax(1) = subplot(4,1,1)
plot(time,IPL,time, IPL_SIGN*1E5, time, IPL_REF_OLD,'LineWidth',2);
title(strcat('Current quench and plateau phase detection. #',num2str(shot)))
grid on;xlim([t_i t_f]);
ax(2) = subplot(4,1,2)
plot(time,DIPL,time,dim*cq_dIpl_threshold,'LineWidth',2);
ylabel('[A/s]');      grid on;xlim([t_i t_f]);
ax(3) = subplot(4,1,3)
plot(time,CQ_DETECT_OUTPUT,'LineWidth',2);
xlabel('time [s]');grid on; xlim([t_i t_f]);
ax(4) = subplot(4,1,4)
plot(time,FBSC3,time,FBSC4,time,FBSC5,time,FBSC6,time,FBSC7,time,FBSC8,'LineWidth',2);
legend('FBSC3','FBSC4','FBSC5','FBSC6','FBSC7','FBSC8','FontSize',14);
grid on; xlim([t_i t_f]);


% ax(4) = subplot(4,1,4)
% plot(time,dim*TRIGGER_threshold_ipl_vloop_activation,time,IPL_ERROR);
% xlabel('time [s]');grid on; xlim([t_i t_f]);


linkaxes(ax,'x');


figure(2)
ax2(1) = subplot(5,2,1)
plot(time,IPL,time, IREF_OUTPUT, time, IPL_REF_OLD,'LineWidth',2);
title(strcat('Ip tracking and Vloop threshold. #',num2str(shot)),'FontSize',14);
legend('ipl','new ire','ipl ref','FontSize',14);
grid on; xlim([t_i t_f]); ylabel('[A]');
ax2(2) = subplot(5,2,3)
plot(time,DIPL,time,dim*cq_dIpl_threshold,'LineWidth',2);
ylabel('[A/s]');      grid on;xlim([t_i t_f]);
legend('dIPL','threshold dIPL','FontSize',14);
ax2(3) = subplot(5,2,5)
plot(time,VLOOP,time,dim*TRIGGER_threshold_vloop*signIPl,'LineWidth',2);
legend('vloop','threshold vloop','FontSize',14);
ylabel('[V]');   grid on;xlim([t_i t_f]);
ax2(4) = subplot(5,2,7)
plot(time,IP_DELTA_OUT,'LineWidth',2);
legend('IP_{DELTA}','FontSize',14);   grid on;xlim([t_i t_f]);
ylabel('[A]');
ax2(5) = subplot(5,2,9)
plot(time, ENABLE_CURRENT_TRACKING,time,ENABLE_DELTA_REF,time,CQ_DETECT_OUTPUT,'LineWidth',2);
legend('ENABLE TRACKING ','ENABLE DELTA REF( enable:5)','CQ DETECT','FontSize',14);   grid on;xlim([t_i t_f]);

xlabel('time [s]');


ax2(6) = subplot(5,2,2)
plot(time,IPL, time, IP_DELTA_PWM, time, IPL_REF_OLD,'LineWidth',2);
title(strcat('Electrical field threshold and square wavePWM trakking.  #',num2str(shot)));
legend('ipl','new ire pwm','ipl ref','FontSize',14);
grid on; xlim([t_i t_f]);
ax2(7) = subplot(5,2,4)
plot(time,DIPL,time,dim*cq_dIpl_threshold,'LineWidth',2);
ylabel('[A/s]');      grid on;xlim([t_i t_f]);
legend('dIPL','threshold dIPL','FontSize',14);
ax2(8) = subplot(5,2,6)
plot(time,VLOOP,time,dim*TRIGGER_threshold_vloop*signIPl,'LineWidth',2);
legend('vloop','threshold vloop','FontSize',14);
ylabel('[V]');   grid on;xlim([t_i t_f]);
ax2(9) = subplot(5,2,8)
plot(time,IREF_LINEAR,'LineWidth',2);
legend('PWM','FontSize',14);   grid on;xlim([t_i t_f]);
ylabel('[A]');
ax2(10) = subplot(5,2,10)
plot(time, ENABLE_CURRENT_TRACKING,time,ENABLE_PWM,time,CQ_DETECT_OUTPUT,'LineWidth',2);
legend('ENABLE TRACKING ','ENABLE PWM','CQ DETECT','FontSize',14);   grid on;xlim([t_i t_f]);
linkaxes(ax2,'x');
xlabel('time [s]');


   
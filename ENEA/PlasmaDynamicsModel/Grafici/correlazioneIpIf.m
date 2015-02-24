clc;clear all;close all;
typeModel = 2;
savePlot=false;

Ts = 0.0005; 
radiusNum=32;
radiusIntPosition=1.22;
ts=0;
ti=0;
%%
Sparo=38864;   
temp1 = CleanFTUdataLOCAL_v4(ti,.8,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp1,Sparo);
Sparo=38865;   
temp2 = CleanFTUdataLOCAL_v4(ti,.8,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp2,Sparo);
Sparo=38869;   
temp3 = CleanFTUdataLOCAL_v4(ti,1.5,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);
Sparo=38519;   
temp4 = CleanFTUdataLOCAL_v4(ti,0.3,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);
Sparo=38513;   
temp5 = CleanFTUdataLOCAL_v4(ti,0.3,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);
Sparo=38863;   
temp6 = CleanFTUdataLOCAL_v4(ti,1.7,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);

Sparo=35965;   
temp11 = CleanFTUdataLOCAL_v4(0.7,1.1,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);

Sparo=36569;   
temp21 = CleanFTUdataLOCAL_v4(0.7,1.3,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);

Sparo=36574;   
temp31 = CleanFTUdataLOCAL_v4(0.7,1.1,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);

Sparo=36634;   
temp41 = CleanFTUdataLOCAL_v4(0.7,1.1,Ts,Sparo,radiusNum,radiusIntPosition);
%graficoIpIfDezDep(temp3,Sparo);

%% rapporto If/Ip
a1=figure('Name','If/Ip');
temp3.mytime=temp3.mytime-0.7;
subplot(2,1,1)
plot(temp1.mytime,temp1.IPLmis, temp2.mytime, temp2.IPLmis, temp3.mytime, temp3.IPLmis,'LineWidth',2); grid on;
grid on; legend('Ip 38864','Ip 38865','Ip 38869','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Plasma current');xlim([0 0.8]);
subplot(2,1,2)
plot(temp1.mytime,abs(temp1.Ifm./temp1.IPLmis), temp2.mytime, abs(temp2.Ifm./temp2.IPLmis), temp3.mytime, abs(temp3.Ifm./temp3.IPLmis),'LineWidth',2); grid on;
grid on; legend('Ip 38864','Ip 38865','Ip 38869','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylim([-0.05 0.15]);
xlim([0 0.8]);
ylabel('','interpreter','latex','fontsize',12);
title('|I_F/I_P|');
saveas(a1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a1'), 'epsc')
saveas(a1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a1'), 'jpeg')
saveas(a1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a1'), 'fig')


a2=figure('Name','If/Ip');
subplot(3,1,1)
plot(temp1.mytime,temp1.Ifm,  temp2.mytime, temp2.Ifm,...
     temp3.mytime, temp3.Ifm,...
     temp1.mytime,temp1.Fcalc,  temp2.mytime, temp2.Fcalc,...
     temp3.mytime, temp3.Fcalc,'LineWidth',2); grid on;
grid on; legend('Ifm 38864','Ifm 38865','Ifm 38869',...
                'Fcalc 38864','Fcalc 38865','Fcalc 38869','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Ifu, Ifcal');
xlim([0 0.8]);
subplot(3,1,2)
plot(temp1.mytime,temp1.dep, temp2.mytime, temp2.dep, temp3.mytime, temp3.dep,'LineWidth',2); grid on;
grid on; legend('dep 38864','dep 38865','dep 38869','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]);
ylabel('','interpreter','latex','fontsize',12);
title('DEP');xlim([0 0.8]);
subplot(3,1,3)
plot(temp1.mytime,temp1.dez, temp2.mytime, temp2.dez, temp3.mytime, temp3.dez,'LineWidth',2); grid on;
grid on; legend('dez 38864','dez 38865','dez 38869','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]);
ylabel('','interpreter','latex','fontsize',12);
title('DEZ');xlim([0 0.8]);
saveas(a2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a2'), 'epsc')
saveas(a2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a2'), 'jpeg')
saveas(a2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\a2'), 'fig')



%%
%% rapporto If/Ip
b1=figure('Name','If/Ip old');
subplot(2,1,1)
temp11.mytime=temp11.mytime+0.2;
plot(temp11.mytime,temp11.IPLmis, temp21.mytime, temp21.IPLmis,...
     temp31.mytime, temp31.IPLmis,temp41.mytime, temp41.IPLmis,'LineWidth',2); grid on;
grid on; legend('Ip 35965','Ip 36569','Ip 36574','Ip 36634','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Plasma current');
xlim([0.9 1.1])
subplot(2,1,2)
plot(temp11.mytime,abs(temp11.Ifm./temp11.IPLmis), temp21.mytime, abs(temp21.Ifm./temp21.IPLmis),...
     temp31.mytime,abs(temp31.Ifm./temp31.IPLmis), temp41.mytime, abs(temp41.Ifm./temp41.IPLmis),'LineWidth',2); grid on;
grid on; legend('35965','36569','36574','36634','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylim([-0.01 0.3]);
xlim([0.9 1.1])
ylabel('','interpreter','latex','fontsize',12);
title('|I_F/I_P|');
saveas(b1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b1'), 'epsc')
saveas(b1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b1'), 'jpeg')
saveas(b1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b1'), 'fig')

b2=figure('Name','If/Ip old');
subplot(3,1,1)
temp11.mytime=temp11.mytime+0.2;
plot(temp11.mytime,temp11.Ifm,  temp21.mytime, temp21.Ifm,...
     temp31.mytime, temp31.Ifm, temp41.mytime, temp41.Ifm,...
     temp11.mytime,temp11.Fcalc,'k-.',temp21.mytime, temp21.Fcalc,'g-.',...
     temp31.mytime, temp31.Fcalc,'r-.', temp41.mytime, temp41.Fcalc,'b-.','LineWidth',2); grid on;
 
grid on; legend('Ifm 35965','Ifm 36569','Ifm 36574','Ifm 36634','Fcalc 35965',...
                'Fcalc 36569','Fcalc 36574','Fcalc 36634','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Ifu, Ifcal');
xlim([0.9 1.1])

subplot(3,1,3)
plot(temp11.mytime,temp11.dez,  temp21.mytime, temp21.dez,...
     temp31.mytime, temp31.dez, temp41.mytime, temp41.dez,'LineWidth',2); grid on;
 
grid on; legend('dez 35965','dez 36569','dez 36574','dez 36634','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]); 
title('DEZ');
ylabel('','interpreter','latex','fontsize',12);
xlim([0.9 1.1]);

subplot(3,1,2)
plot(temp11.mytime,temp11.dep,  temp21.mytime, temp21.dep,...
     temp31.mytime, temp31.dep, temp41.mytime, temp41.dep,'LineWidth',2); grid on; 
grid on; legend('dep 35965','dep 36569','dep 36574','dep 36634','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]);
title('DEP')
ylabel('','interpreter','latex','fontsize',12);
xlim([0.9 1.1]);

saveas(b2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b2'), 'epsc')
saveas(b2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b2'), 'jpeg')
saveas(b2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\b2'), 'fig')







%% rapporto If/Ip
c1=figure('Name','If/Ip bis');
subplot(2,1,1)
plot(temp4.mytime,temp4.IPLmis, temp5.mytime, temp5.IPLmis,'LineWidth',2); grid on;
grid on; legend('Ip 38519','Ip 38513','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Plasma current');xlim([0 0.3])
subplot(2,1,2)
plot(temp4.mytime,abs(temp4.Ifm./temp4.IPLmis), temp5.mytime, abs(temp5.Ifm./temp5.IPLmis), 'LineWidth',2); grid on;
grid on; legend('38519','38513','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylim([-0.04 0.5]);
xlim([0 0.3])
ylabel('','interpreter','latex','fontsize',12);
title('|I_F/I_P|');

saveas(c1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c1'), 'epsc')
saveas(c1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c1'), 'jpeg')
saveas(c1,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c1'), 'fig')
%%
c2=figure('Name','If dez dep');
subplot(3,1,1)
plot(temp4.mytime,temp4.Ifm,   temp5.mytime, temp5.Ifm,...
     temp4.mytime,temp4.Fcalc,'k-.', temp5.mytime, temp5.Fcalc,'r-.', 'LineWidth',2); grid on;
grid on; legend('Ifu 38519','Ifu 38513','Ifcal 38519','Ifcal 38513','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Ifu, Ifcal');
xlim([0 0.3])
subplot(3,1,2)
plot(temp4.mytime,temp4.dep, temp5.mytime, temp5.dep,'LineWidth',2); grid on;
grid on; legend('dep 38519','dep 38519','dez 38519','dez 38519','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]);
xlim([0 0.3])
ylabel('Web','interpreter','latex','fontsize',12); title('DEP');
subplot(3,1,3)
plot(temp4.mytime,temp4.dez, temp5.mytime, temp5.dez,'LineWidth',2); grid on;
grid on; legend('dez 38519','dez 38519','Location','west');
xlabel('time[s]','interpreter','latex','fontsize',12);
%ylim([-0.1 0.1]);
xlim([0 0.3])
ylabel('','interpreter','latex','fontsize',12);
title('DEZ');

saveas(c2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c2'), 'epsc')
saveas(c2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c2'), 'jpeg')
saveas(c2,strcat('C:\Users\mateusz\Google Drive\Dottorato2014\ENEA\plasmaElongation\figs\c2'), 'fig')
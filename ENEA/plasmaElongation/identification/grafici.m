%% Identificazione dez
%% 6/6/1014

clc;clear all;close all;
Ts = 0.0005; 
radiusNum=32;
radiusIntPosition=1.22;


%37862 ... %FTU-X
%37864 ... %FTU-X
%37869 ... %FTU-X
%37926 ... %FTU-X
%37929 ... %FTU-X
%37930 ... %FTU-X nuovi guadagni PIDH
%37931 ... %FTU-X nuovi guadagni PIDH
%37932 ... %FTU-X nuovi guadagni PIDH, instabilità verticale, cambiati i riferimenti, messi più esterni e quindi perso il plasma verticalmente
%38298 ... %oscillazione + disruzione vertical FTU-x
%37869 ... %FTU-X perdita rampa down
%37870 ... %FTU-X oscillazioni riprese
    


ds=0.3;
shot=37870;    
ts=0
ti37870=ts;
temp37870 = CleanFTUdataLOCAL_v4(ti37870,1.352+ds,Ts,shot,radiusNum,radiusIntPosition);

shot=37869;
ti37869=ts;
temp37869=  CleanFTUdataLOCAL_v4(ti37869,1.48+ds,Ts,shot,radiusNum,radiusIntPosition);

shot=38298; 
ti38298=ts;
temp38298 =  CleanFTUdataLOCAL_v4(ti38298,0.9+ds,Ts,shot,radiusNum,radiusIntPosition);

shot=37932; 
ti37932=ts;
temp37932 =  CleanFTUdataLOCAL_v4(ti37932,0.885+ds,Ts,shot,radiusNum,radiusIntPosition);


%% Grafici 
coilH_IHcal=figure('Name','Coil H');
subplot(2,2,1);
plot(temp37870.mytime,[temp37870.Hcalc'],'k','LineWidth',1); 
grid on; legend('I_H'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37870)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(2,2,2);
plot(temp37869.mytime,[temp37869.Hcalc'],'k','LineWidth',1); 
grid on; legend('I_H'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37869)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(2,2,3);
plot(temp38298.mytime,[temp38298.Hcalc'],'k','LineWidth',1); 
grid on; legend('I_H'); title(strcat('Corrente H calcolata. Sparo: ', num2str(38298)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(2,2,4);
plot(temp37932.mytime,[temp37932.Hcalc'],'k','LineWidth',1); 
grid on; legend('I_H'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37932)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
saveas(coilH_IHcal,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_IHcal'), 'fig')
saveas(coilH_IHcal,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_IHcal'), 'epsc')
saveas(coilH_IHcal,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_IHcal'), 'jpeg')


coilH_Dez=figure('Name','Dez');
subplot(2,2,1);
plot(temp37870.mytime,[temp37870.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37870)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(2,2,2);
plot(temp37869.mytime,[temp37869.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37869)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(2,2,3);
plot(temp38298.mytime,[temp38298.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(38298)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(2,2,4);
plot(temp37932.mytime,[temp37932.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37932)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_Dez'), 'fig')
saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_Dez'), 'epsc')
saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilH_Dez'), 'jpeg')


coilH_IP=figure('Name','IP')
subplot(2,2,1);
plot(temp37870.mytime,[temp37870.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corrente di plasma. Sparo: ', num2str(37870)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,2);
plot(temp37869.mytime,[temp37869.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corrente di plasma. Sparo: ', num2str(37869)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,3);
plot(temp38298.mytime,[temp38298.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corrente di plasma. Sparo: ', num2str(38298)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,4);
plot(temp37932.mytime,[temp37932.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corrente di plasma. Sparo: ', num2str(37932)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'fig')
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'epsc')
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'jpeg')



%%
sparo38298=figure('Name','IP38298')
subplot(3,1,1);
plot(temp38298.mytime,[temp38298.Hcalc' temp38298.Hmis'],'LineWidth',1); 
grid on; legend('I_HCal','I_HMis'); title(strcat('Corrente H calcolata. Sparo: ', num2str(38298)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(3,1,2);
plot(temp38298.mytime,[temp38298.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(38298)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(3,1,3);
plot(temp38298.mytime,[temp38298.Ip' temp38298.vloop'*1e4],'LineWidth',1); 
grid on; legend('Ip','V_{loop}1e4'); title(strcat('Corrente di plasma. Sparo: ', num2str(38298)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(sparo38298,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo38298'), 'fig')
saveas(sparo38298,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo38298'), 'epsc')
saveas(sparo38298,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo38298'), 'jpeg')
%%
%%
sparo37932=figure('Name','IP37932')
subplot(3,1,1);
plot(temp37932.mytime,[temp37932.Hcalc' temp37932.Hmis'],'LineWidth',1); 
grid on; legend('I_HCal','I_HMis'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37932)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(3,1,2);
plot(temp37932.mytime,[temp37932.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37932)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(3,1,3);
plot(temp37932.mytime,[temp37932.Ip' temp37932.vloop'*1e4],'LineWidth',1); 
grid on; legend('Ip','V_{loop}1e4'); title(strcat('Corrente di plasma. Sparo: ', num2str(37932)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(sparo37932,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37932'), 'fig')
saveas(sparo37932,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37932'), 'epsc')
saveas(sparo37932,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37932'), 'jpeg')
%%
sparo37869=figure('Name','IP37869')
subplot(3,1,1);
plot(temp37869.mytime,[temp37869.Hcalc' temp37869.Hmis'],'LineWidth',1); 
grid on; legend('I_HCal','I_HMis'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37869)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(3,1,2);
plot(temp37869.mytime,[temp37869.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37869)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(3,1,3);
plot(temp37869.mytime,[temp37869.Ip' temp37869.vloop'*1e4],'LineWidth',1); 
grid on; legend('Ip','V_{loop}1e4'); title(strcat('Corrente di plasma. Sparo: ', num2str(37869)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(sparo37869,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37869'), 'fig')
saveas(sparo37869,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37869'), 'epsc')
saveas(sparo37869,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37869'), 'jpeg')
%%
sparo37870=figure('Name','IP37870')
subplot(3,1,1);
plot(temp37870.mytime,[temp37870.Hcalc' temp37870.Hmis'],'LineWidth',1); 
grid on; legend('I_HCal','I_HMis'); title(strcat('Corrente H calcolata. Sparo: ', num2str(37870)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[A]','interpreter','latex','fontsize',10);
subplot(3,1,2);
plot(temp37870.mytime,[temp37870.dez'],'k','LineWidth',1); 
grid on; legend('dez'); title(strcat('Flusso verticale. Sparo: ', num2str(37870)));
xlabel('tempo[s]','interpreter','latex','fontsize',10);
ylabel('[Web]','interpreter','latex','fontsize',10);
subplot(3,1,3);
plot(temp37870.mytime,[temp37870.Ip' temp37870.vloop'*1e4],'LineWidth',1); 
grid on; legend('Ip','V_{loop}1e4'); title(strcat('Corrente di plasma. Sparo: ', num2str(37870)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(sparo37870,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37870'), 'fig')
saveas(sparo37870,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37870'), 'epsc')
saveas(sparo37870,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\sparo37870'), 'jpeg')
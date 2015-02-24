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
grid on; legend('Ip'); title(strcat('Corente di plasma. Sparo: ', num2str(37870)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,2);
plot(temp37869.mytime,[temp37869.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corente di plasma. Sparo: ', num2str(37869)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,3);
plot(temp38298.mytime,[temp38298.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corente di plasma. Sparo: ', num2str(38298)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
subplot(2,2,4);
plot(temp37932.mytime,[temp37932.Ip'],'k','LineWidth',1); 
grid on; legend('Ip'); title(strcat('Corente di plasma. Sparo: ', num2str(37932)));
ylabel('[A]','interpreter','latex','fontsize',10);
xlabel('tempo[s]','interpreter','latex','fontsize',10);
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'fig')
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'epsc')
saveas(coilH_IP,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\Ip'), 'jpeg')

order=[1 1 3];                      % [Ny Nu Nx]
t=0;
break
%%

init=zeros(3,1);
initCell = {init(1,:) init(2,:) init(3,:)};
initialStates = struct('Name',  {'x1' 'x2' 'x3'},...
                       'Unit',  {'' '' ''},...
                       'Value',  initCell,...
                       'Minimum', {-Inf -Inf -Inf},...
                       'Maximum',{Inf Inf Inf},...
                       'Fixed',{0 0 0});
par=[1.339e2 1.522e4 5.512e5 0 0 -2.753E2]; % da definire

% 
par=[8.3331e+000
     1.9691e+002
    -3.1590e+000
    -1.3953e-001
     2.7413e+001
     1.3575e-001];


 par=[2.1452e2
     3.7642e3
     1.6938e4
    -1.3884e-1
     -1.5104
      -3.8962e-1];

model = idnlgrey('modelCoilH',order,par);%,initialStates,t,'Name', 'Model Coil H Identification');

%% parameter
setpar(model, 'Fixed', {0 0 0 0 0 0}); 
setpar(model, 'Fixed', {1 1 1 1 1 1}); 
%setpar(model, 'Minimum', {0 1 20 0 -Inf 0 0}); 
%setpar(model, 'Maximum', {Inf 40 100 1 Inf}); 
z1 = temp37870.zIddData;z1Sim = temp37870.zIddDataSim;
z2 = temp37869.zIddData;z2Sim = temp37869.zIddDataSim;
z3 = temp38298.zIddData;z3Sim = temp38298.zIddDataSim;
z4 = temp37932.zIddData;z4Sim = temp37932.zIddDataSim;

%% parameter
model.Algorithm.SimulationOptions.Solver= 'ode4';
model.Algorithm.Display='On';
%model.Algorithm.Tolerance=1e-7;
%model.Algorithm.MaxIter=20;
%%
display('Merge data ');
z = merge(z1,z2,z3,z4); %multiexperiment   

%identificazione del modello
display('Start Identification');
modelId = pem(z, model);
param = getpar(modelId,'Value')
%%
[YH, FIT, X0]=compare(modelId,z);
%%
close all;
y=get(YH{1},'OutputData')
coilH_Dez=figure('Name',strcat('Dez shot'))
sparo=37870;
subplot(2,2,1)
%coilH_Dez=figure('Name',strcat('Dez shot',num2str(sparo)))
plot(y{1},'k','LineWidth',1.5); hold on; plot(get(z1,'OutputData'),'--r','LineWidth',0.5)
grid on; legend('Dez identificato','Dez misurato','Location','Best'); 
title(strcat('Identificazione del flusso verticale. Sparo: ', num2str(sparo)));
ylabel('[Wb]','interpreter','latex','fontsize',10);

% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'fig')
% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'epsc')
% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'jpeg')

sparo=37869;
%coilH_Dez=figure('Name',strcat('Dez shot',num2str(sparo)))
subplot(2,2,2)
plot(y{2},'k','LineWidth',1.5); hold on; plot(get(z2,'OutputData'),'--r','LineWidth',0.5)
grid on; legend('Dez identificato','Dez misurato','Location','Best'); 
title(strcat('Identificazione del flusso verticale. Sparo: ', num2str(sparo)));
ylabel('[Wb]','interpreter','latex','fontsize',10);
% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'fig')
% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'epsc')
% saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'jpeg')

sparo=38298;
%coilH_Dez=figure('Name',strcat('Dez shot',num2str(sparo)))
subplot(2,2,3)
plot(y{3},'k','LineWidth',1.5); hold on; plot(get(z3,'OutputData'),'--r','LineWidth',0.5)
grid on; legend('Dez identificato','Dez misurato','Location','Best'); 
title(strcat('Identificazione del flusso verticale. Sparo: ', num2str(sparo)));
ylabel('[Wb]','interpreter','latex','fontsize',10);
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'fig')
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'epsc')
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'jpeg')

sparo=37932;
%coilH_Dez=figure('Name',strcat('Dez shot',num2str(sparo)))
subplot(2,2,4)
plot(y{4},'k','LineWidth',1.5); hold on; plot(get(z4,'OutputData'),'--r','LineWidth',0.5)
grid on; legend('Dez identificato','Dez misurato','Location','Best'); 
title(strcat('Identificazione del flusso verticale. Sparo: ', num2str(sparo)));
ylabel('[Wb]','interpreter','latex','fontsize',10);
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'fig')
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'epsc')
%saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo',num2str(sparo)), 'jpeg')

saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo'), 'fig')
saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo'), 'epsc')
saveas(coilH_Dez,strcat('D:\Dropbox\TesiRunaway\Latex\imgs\coilH\coilHSparo'), 'jpeg')


%% modello
% 
clc;
param = getpar(modelId,'Value');    
p= cell2mat(param);
A=[-p(1)       1        0;
   -p(2)       0        1;
   -p(3)       0        0]

B=[p(4); p(5); p(6)]
autovalori = eig(A)
C=[1 0 0 ]
D=[0]
coilHSystem = ss(A,B,C,D);
[Ad,Bd,Cd,Dd] = ssdata(c2d(coilHSystem,1e-3,'zoh'));

%     [ 2.1452e+002]
%     [ 3.7642e+003]
%     [ 1.6938e+004]
%     [-1.3884e-001]
%     [-1.5104e+000]
%     [-3.8962e-001]

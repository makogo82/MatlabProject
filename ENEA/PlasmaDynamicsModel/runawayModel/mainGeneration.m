clear all; clc; close all;


data = load('data37714.mat');  %dati zana
a1=data.a;
b1=data.b;
c1=data.c;
d1=data.d;

b=ones(1,20)/20;
t = (data.time);
n = data.n(2:size(t));
TeMax = data.TeMax(2:size(t)); %- 'TeMax' (central electron temperature) is in keV
neMax = data.neMax(2:size(t));
Ell = data.Ell(2:size(t));
Spri = data.Spri(2:size(t));
Zeff = data.Zeff(2:size(t));

title('$\int_{}^{eee}$','interpreter','latex')

SpriCalc = Dreicer( neMax,Zeff,TeMax,Ell);
plot(t(2:size(t)),SpriCalc.Spri*.4e-12,t(2:size(t)),Spri,'r');
ylim([0 1e+8]); grid on;

%%
close all;
figure(1)
subplot(2,1,1)
plot(t(2:size(t)),Spri); hold on;
subplot(2,1,2)
plot(t(2:size(t)),SpriCalc.Spri);

break;
%%
dati = load('shot_37714.mat'); %dati cluster
NEU213y = dati.NEU213.y;
NEU213x = dati.NEU213.x;
BF3x = dati.BF3.x;
BF3y = dati.BF3.y;
IPL = dati.IPLmis.y;
time = dati.IPLmis.x
time = dati.IPLmis.x;
Vloop = filtfilt(b,1,dati.Vloop.y);
zeff = dati.ZEFF.y;
time = dati.ZEFF.x;
close all
plot(time,zeff)

DensitySINVEA = dati.DensitySINVEA;
d=[]

for t = 1:61
 d(t) = max(DensitySINVEA.z(:,t));
end
time =(data.time)
figure('Name','de')
plot(DensitySINVEA.y,d,time(2:size(time)),neMax,'r')

figure('Name','Te')
plot(dati.TeMax.x,dati.TeMax.y,'k'); hold on; plot(time(2:size(time)),TeMax);

mytime=[DensitySINVEA.y(1):0.005:DensitySINVEA.y(length(DensitySINVEA.y))];
densitaMax = interp1(DensitySINVEA.y,d,mytime,'spline');
figure('Name','MyDensity')
plot(mytime,densitaMax)

SpriCalc = Dreicer( neMax,Zeff,TeMax,Ell);


%ZEFF
%TeMax
%DensitySINVEA
%DensitySIDE3



%%
SpriMy=Dreicer( neMax,dati.ZEFF.y,TeMax,dati.VloopFilter.y);


break;

a=figure('Name','Zeff')
plot(t,dati); title('Zeff');

NEU213yFilter = filtfilt(b,1,NEU213y); 
BF3yFilter = filtfilt(b,1,BF3y); 

NEU213AndBF3 = figure('Name','NEU213AndBF3');
%subplot(2,1,1)
semilogy(NEU213x,NEU213yFilter,'k',BF3x,BF3yFilter,'g','LineWidth',2); grid on; title('log_1_0');
grid on;
xlabel('time(s)');
legend('NEU213','BF3','Location','northwest'); 
xlim([0,1.5]); plotParam;
ylim([3e9 1e13]);
% subplot(2,1,2)
% plot(NEU213x,NEU213yFilter,'k',BF3x,BF3y,'g','LineWidth',2); grid on; title('normal');
% grid on
% legend('NEU213y','BF3','Location','northwest'); plotParam;
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'epsc')
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'jpeg')
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'fig')

%%


%% controllo se ci sono dei Nan
for j=1:size(Spri)
    if (isnan(Spri(j)))
        Spri(j)=Spri(j+1);
    end
end
t = t(2:size(t));    
%% interpolazione dei segnali 
newTime = [0.35:0.0005:1.4]';
nN = interp1(t,n,newTime);
TeMaxN = interp1(t,TeMax,newTime);
neMaxN = interp1(t,neMax,newTime);
EllN = interp1(t,Ell,newTime);
SpriN = interp1(t,Spri,newTime);
iplN = interp1(time,IPL,newTime);
vloopN = interp1 ( time , Vloop, newTime);
NEU213yN= interp1(NEU213x,NEU213yFilter,newTime);
BF3yN =   interp1(BF3x, BF3y, newTime);





%% Integrale
for i = 2: size(SpriN,1)
 IntSpri(i) = trapz(newTime(1:i),SpriN(1:i));
end
%%
figure(1);
nsubplot = 3;
subplot(nsubplot,1,1)
plot(newTime,TeMaxN,'LineWidth',2); grid on; title('TeMax');ylabel('[keV]');
xlim([newTime(1),newTime(size(newTime,1)-1)]); plotParam;
subplot(nsubplot,1,2)
plot(newTime,EllN,'LineWidth',2); grid on;  title('E_{||}'); ylabel('[V/m]');
xlim([newTime(1),newTime(size(newTime,1)-1)]); plotParam;
subplot(nsubplot,1,3)
plot(newTime,neMaxN,'LineWidth',2); grid on; title('density');
xlabel('time(s)');ylabel('[m^{-3}]');
xlim([newTime(1),newTime(size(newTime,1)-1)]); plotParam;
%%
figure('Name','SpriN and IntSpri')
ax2(1)= subplot(2,1,1)
plot(newTime,SpriN,'LineWidth',2); grid on; title('S_{pri}'); plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ax2(2)= subplot(2,1,2)
plot(newTime,IntSpri,'LineWidth',2); grid on; title('S_{int}');
linkaxes(ax2,'x');
xlim([newTime(1),newTime(size(newTime,1)-1)]); plotParam;
%%
figure(3)
subplot(2,1,1)
plot(newTime,iplN,'LineWidth',2); grid on; title('IPL'); plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
subplot(2,1,2)
plot(newTime,vloopN,'LineWidth',2); grid on; title('vloopN'); plotParam;
figure(4)
plot(newTime,NEU213yN,newTime,NEU213yN-BF3yN,newTime,IntSpri*6e4,'LineWidth',2); plotParam;
grid on; title('NEU213-BF3'); legend('NEU213','NEU213-BF3','Int(S_{pri})*6E4','Location','northwest');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
%%
NEU213 = figure('Name','Diff')
ax(1)=subplot(2,1,1)
plot(newTime,NEU213yN-BF3yN,newTime,IntSpri*6e4,'LineWidth',2); plotParam;
grid on; title('NEU213-BF3'); legend('NEU213-BF3','Int(S_{pri})*6E4');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ylim([-1e10,2e11]);
ylabel('[a.u]');
ax(2)=subplot(2,1,2)
plot(newTime,mypseudo_derivative(newTime,NEU213yN-BF3yN,40,150),newTime,mypseudo_derivative(newTime,IntSpri*6e4,40,150),'LineWidth',2);
grid on; title('d(NEU213-BF3)/dt and S_{pri}'); 
legend('d(NEU213-BF3)/dt','S_{pri}','Location','northwest');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
xlabel('time(s)');
ylabel('[a.u]');
linkaxes(ax,'x');
plotParam;
saveas(NEU213,strcat('plot\NEU213'), 'epsc')
saveas(NEU213,strcat('plot\NEU213'), 'jpeg')
saveas(NEU213,strcat('plot\NEU213'), 'fig')


%%
p0=130
p11=0.275/0.1843;
figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,neMaxN*2E-21,newTime,EllN,'LineWidth',2); grid on; title('Ell, n_e');
legend('n_e*1E-8','Ell*1E13');plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ax2(2)= subplot(3,1,2)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN),10,60),newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*neMaxN*2E-21+EllN),'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; title('d(NEU213-BF3)/dt'); legend('NEU213-BF3');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
xlabel('time(s)','interpreter','latex')
ax2(3)= subplot(3,1,3)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; title('d(NEU213-BF3)/dt'); legend('NEU213-BF3');

linkaxes(ax2,'x');
xlim([0.2 0.6])
plotParam;


%% sistema MISO
% DAT = iddata(Y,U,Ts)
z = iddata(NEU213yN-BF3yN,[neMaxN EllN],0.0005);
z.OutputName = 'n_e' ;
z.OutputUnit = 'count' ;
z.InputName = {'Ell' 'Electron Density'};
z.InputUnit = {'V/m' 'm-3'};


%primo tentativo n4sid
%sysN4sid = n4sid(z,1);
%figure('Name','n4sid');
%compare(z,sysN4sid);



%%
p0=400
p11=1.00*0.2411/0.1941;
figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,neMaxN*6E-21,newTime,EllN,'LineWidth',2); grid on; title('Ell, n_e');
legend('n_e*1E-8','Ell*1E13');plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ax2(2)= subplot(3,1,2)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN),10,60),newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*neMaxN*6E-21+EllN),'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; title('d(NEU213-BF3)/dt'); legend('NEU213-BF3');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
xlabel('time(s)','interpreter','latex')
ax2(3)= subplot(3,1,3)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; title('d(NEU213-BF3)/dt'); legend('NEU213-BF3');

linkaxes(ax2,'x');
xlim([0.35 0.75])
plotParam;



break

%% PEM

order=[1 2 1]; % [Ny Nu Nx]
t=0;
initialStates = struct('Name',  {'x1'},...
                       'Unit',  {''},...
                       'Value',  {IntSpri(1)},...
                       'Minimum', {-Inf},...
                       'Maximum',{Inf},...
                       'Fixed',{1});
                               
                   
%par=[-4.8635e-10 0 9.9219]
par=[7.1767e-08 9.9998e+12 4.7989]
%par=[7.0703e-08 9.9998e+12 4.4088]
model = idnlgrey('runaway',order,par);%,initialStates,t,'Name', 'Model Coil H Identification');
setpar(model, 'Fixed', {1 1 1}); 
setpar(model, 'Minimum', {0 0 0});
%setpar(model, 'Maximum', {Inf Inf 0.001});

%parameter
model.Algorithm.SimulationOptions.Solver= 'ode45';
model.Algorithm.Display='On';
model.Algorithm.Tolerance=1e-8;
model.Algorithm.MaxIter=40;
%
display('Merge data ');
z = merge(z); %multiexperiment   
%identificazione del modello
display('Start identification grey-box model');
modelId = pem(z, model);

figure('Name','pem');
%%
[YH, FIT, X0] = compare(z,modelId);
param = getpar(modelId,'Value')
y=get(YH,'OutputData');
figure('Name','compareData')
%coilH_Dez=figure('Name',strcat('Dez shot',num2str(sparo)))
plot(newTime,y,'b',newTime,NEU213yN-BF3yN,'k',newTime,IntSpri*6e4,'r','LineWidth',2); plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
grid on; legend('model (PEM)','Experimental data (NEU213-BF3)','Dreicer Eq. (IntSpri*6e4)','Location','northwest'); 
title(strcat('Runaway generation model Identification. Fitting: ', num2str(FIT),'%'));
xlabel('time(s)','interpreter','latex','fontsize',10);



NEU213 = figure('Name','Diff')
ax(1)=subplot(2,1,1)
plot(newTime,EllN*9.9E12,newTime,neMaxN*7.17E-8,'LineWidth',2);
grid on; title('Input'); legend('Ell*9.9E12','n_e*7.17E-8','Location','northwest');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
plotParam;
ax(2)=subplot(2,1,2)
plot(newTime,y,'b',newTime,(NEU213yN-BF3yN),'k',newTime,IntSpri*6e4,'r','LineWidth',1.5);
grid on; title('Output'); legend('model (PEM)','(NEU213-BF3)','IntSpri*6E4','Location','northwest'); 
xlim([newTime(1),newTime(size(newTime,1)-1)]);
plotParam;
xlabel('$time(s)$','interpreter','latex')
ylabel('$n_{RE}$[count]','interpreter','latex')
linkaxes(ax2,'x');



%%
p0=90
p11=1.00*0.2411/0.1941;
%p0=400
%p11=1.00*0.2411/0.1941;
figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,neMaxN*6E-21,newTime,EllN,'LineWidth',2); grid on; title('Ell, n_e');
legend('n_e*1E-8','Ell*1E13');plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ax2(2)= subplot(3,1,2)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN),10,40),newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*neMaxN*6E-21+EllN),'g',newTime,SpriN*9.5E4,'r','LineWidth',2);
grid on; title(''); 
legend('d(NEU213-BF3)/dt','model', 'SpriN*9.5E4');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
xlabel('time(s)','interpreter','latex')
ax2(3)= subplot(3,1,3)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,NEU213yN-BF3yN,newTime,IntSpri*8e4,'r','LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; 
legend('NEU213-BF3','IntSpri*9.5e4'); 

linkaxes(ax2,'x');
xlim([0.35 1.4])
plotParam;



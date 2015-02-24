clear all; 
clc; 
close all;



data37405 = load('data37405.mat');  %dati zana
%%
newTime = [0.1:0.0005:0.7]';
b=ones(1,20)/20;
t = (data37405.t);
n = data37405.n(2:size(t));
TeMax = data37405.TeMax(2:size(t)); %- 'TeMax' (central electron temperature) is in keV
neMax = data37405.neMax(2:size(t));
Ell = data37405.Ell(2:size(t));
Spri = data37405.Spri(2:size(t));

shot37405 = load('shot_37405.mat'); %dati cluster
dati = shot37405.Data;

dens =  dati.esidens.y;
denst = dati.esidens.x;

NEU213y = dati.NEU213.y;
NEU213x = dati.NEU213.x;
BF3x = dati.BF3.x;
BF3y = dati.BF3.y;
IPL = shot37405.Data.IPLmis.y;
time = shot37405.Data.IPLmis.x;
Vloop = filtfilt(b,1,shot37405.Data.Vloop.y);


NEU213yFilter = filtfilt(b,1,NEU213y); 
BF3yFilter = filtfilt(b,1,BF3y); 

NEU213AndBF3 = figure('Name','NEU213AndBF3');
semilogy(NEU213x,NEU213yFilter,'k',BF3x,BF3yFilter,'g','LineWidth',2); grid on; title('log_1_0');
grid on; legend('NEU213','BF3','Location','northwest'); plotParam;
xlim([0,1]);ylim([1e9,1e12]); 
xlabel('time(s)');
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'epsc')
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'jpeg')
saveas(NEU213AndBF3,strcat('plot\NEU213AndBF3'), 'fig')



%%


%% controllo se ci sono dei Nan
for j=1:size(Spri)
    SpriNum(j)= str2double(Spri(j));
    if (isnan(SpriNum(j)))
        if(j>1)
            SpriNum(j)=SpriNum(j-1);
        end
        if(j==1)
            SpriNum(j)=0;
        end
    end
end
t = t(2:size(t));    
%% interpolazione dei segnali 

nN = interp1(t,n,newTime); 
TeMaxN = interp1(t,TeMax,newTime);
neMaxN = interp1(t,neMax,newTime); 
EllN = interp1(t,Ell,newTime);
SpriN = interp1(t,SpriNum,newTime);
iplN = interp1(time,IPL,newTime);
vloopN = interp1 ( time , Vloop, newTime);
NEU213yN= interp1(NEU213x,NEU213yFilter,newTime);
BF3yN =   interp1(BF3x, BF3y, newTime);
densN = interp1(denst,dens,newTime);

%% Integrale
for i = 2: size(SpriN,1)
 IntSpri(i) = trapz(newTime(1:i),SpriN(1:i));
end
%%
figure(1);
nsubplot = 4;
subplot(nsubplot,1,1)
plot(newTime,TeMaxN,'LineWidth',2); grid on; title('TeMax');ylabel('[keV]');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
subplot(nsubplot,1,2);plotParam;
plotParam;
plot(newTime,EllN,'LineWidth',2); grid on;  title('E_{||}');ylabel('[V/m]');
xlim([newTime(1),newTime(size(newTime,1)-1)]);plotParam;
subplot(nsubplot,1,3)
plot(newTime,neMaxN,'LineWidth',2); grid on; title('density');ylabel('[m^{-3}]');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
plotParam;
subplot(nsubplot,1,4)
plot(newTime,vloopN,'LineWidth',2); grid on; title('vloop');
xlabel('time(s)','interpreter','latex');ylabel('[V]');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
plotParam;
%%
p0=130
p11=1.1*0.275/0.1843;
figure37405=figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,neMaxN*2E-21,newTime,EllN,'LineWidth',2); grid on;
legend('n_e*2E-21','Ell');plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);plotParam;
ax2(2)= subplot(3,1,2);
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
dS=mypseudo_derivative(newTime,(NEU213yN-BF3yN),5,30);
plot(newTime,dS,newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*neMaxN*2E-21+EllN),'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; legend('d(NEU213-BF3)/dt','model');
xlim([newTime(1),newTime(size(newTime,1)-1)]);plotParam;
ax2(3)= subplot(3,1,3)
countRE = NEU213yN-BF3yN;
x0=countRE(find(newTime > 0.2485 ,1));
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%%
[T,Y] = ode45(@(t,y) model(t,y,newTime,neMaxN, EllN), [0.2485 0.7], x0);
plot(newTime,NEU213yN-BF3yN,T,Y,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; legend('NEU213-BF3','model');
xlabel('time(s)','interpreter','latex')
linkaxes(ax2,'x');
xlim([0.2485 0.7])
plotParam;

saveas(figure37405,strcat('plot\figure37405'), 'jpeg')
saveas(figure37405,strcat('plot\figure37405'), 'fig')
saveas(figure37405,strcat('plot\figure37405'), 'png')
saveas(figure37405,strcat('plot\figure37405'), 'epsc')








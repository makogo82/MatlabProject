clear all; 
clc; 
close all;


%%
newTime = [0.12:0.0005:0.9]';
b=ones(1,20)/20;


%%
dati = shot.Data;
NEU213y = dati.NEU213.y;
NEU213x = dati.NEU213.x;
BF3x = dati.BF3.x;
BF3y = dati.BF3.y;
IPL = shot.Data.IPLmis.y;
time = shot.Data.IPLmis.x;
Vloop = filtfilt(b,1,shot.Data.Vloop.y);
dens =  dati.esidens.y;
denst = dati.esidens.x;

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


%% interpolazione dei segnali 
iplN = interp1(time,IPL,newTime);
vloopN = interp1 ( time , Vloop, newTime);
NEU213yN= interp1(NEU213x,NEU213yFilter,newTime);
BF3yN =   interp1(BF3x, BF3y, newTime);
densN = interp1(denst,dens,newTime);

%%
p0=5
p11=1.5*0.275/0.1843;
p11 =2.405/0.8363;
figure37405=figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,densN*2E-20,newTime,vloopN,'LineWidth',2); grid on;
legend('n_e*2E-20','Ell');plotParam;
xlim([newTime(1),newTime(size(newTime,1)-1)]);plotParam;
ax2(2)= subplot(3,1,2);

plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN),10,50),newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*densN*2E-20+vloopN),'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; legend('d(NEU213-BF3)/dt','model');
xlim([newTime(1),newTime(size(newTime,1)-1)]);plotParam;
ax2(3)= subplot(3,1,3)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');

plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; title('NEU213-BF3'); legend('NEU213-BF3');
xlabel('time(s)','interpreter','latex')
linkaxes(ax2,'x');
xlim([0.1 0.9])
plotParam;

saveas(figure37405,strcat('plot\figure37405'), 'jpeg')
saveas(figure37405,strcat('plot\figure37405'), 'fig')
saveas(figure37405,strcat('plot\figure37405'), 'png')
saveas(figure37405,strcat('plot\figure37405'), 'epsc')

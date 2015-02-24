clear all; 
clc; 
close all;

%%
newTime = [0.7:0.0005:1.2]';
data = load('data38358.mat');  %dati zana
b=ones(1,20)/20;
%1280
t = (data.ts(1:1280));
n = data.ne1019m3(2:size(t));
TeMax = data.TeMaxkeV(2:size(t)); %- 'TeMax' (central electron temperature) is in keV
neMax = data.neMaxm3(2:size(t));
Ell = data.EllVm(2:size(t));
Spri = data.Spriau(2:size(t));

dati = load('shot_38358.mat'); %dati cluster

NEU213y = dati.NEU213.y;
NEU213x = dati.NEU213.x;
BF3x = dati.BF3.x;
BF3y = dati.BF3.y;
IPL = dati.IPLmis.y;
time = dati.IPLmis.x;
Vloop = filtfilt(b,1,dati.Vloop.y);

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
%vloopN = interp1 ( time , Vloop, newTime);
NEU213yN= interp1(NEU213x,NEU213yFilter,newTime);
BF3yN =   interp1(BF3x, BF3y, newTime);

%% Integrale
for i = 2: size(SpriN,1)
 IntSpri(i) = trapz(newTime(1:i),SpriN(1:i));
end

%% generation

%%
p0=90
p11=1.00*0.2411/0.1941;
%p0=400
%p11=1.00*0.2411/0.1941;
figure('Name','Ell-n_e-NEU213-BF3')
ax2(1)= subplot(3,1,1)
plot(newTime,neMaxN,newTime,EllN,'LineWidth',2); grid on; title('Ell, n_e');
legend('n_e*1E-8','Ell*1E13');plotParam;
title('SHOT # 37714')
xlim([newTime(1),newTime(size(newTime,1)-1)]);
ax2(2)= subplot(3,1,2)
%plot(newTime,NEU213yN-BF3yN,'LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
%grid on; title('NEU213-BF3'); legend('NEU213-BF3');
plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN),10,40),newTime,...
     p0*(NEU213yN-BF3yN).*(-p11*neMaxN*6E-21+EllN),'g',newTime,SpriN*9.5E4,'r','LineWidth',2);
grid on; title(''); 
legend('d(NEU213-BF3)/dt','model', 'Spri*9.5E4');
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
%%

%%
p0=90
p11=1.00*0.2411/0.1941;
%p0=400
%p11=1.00*0.2411/0.1941;
figure37714=figure('Name','A')
ax2(2)= subplot(2,1,1)
plot(newTime,mypseudo_derivative(newTime,(NEU213yN-BF3yN)*5e3,10,40),newTime,SpriN,'r','LineWidth',2);
title('SHOT # 37714'); grid on;
legend('d(NEU213-BF3)/dt', 'Spri','interpreter','latex','Location','northwest');
xlim([newTime(1),newTime(size(newTime,1)-1)]);
plotParam;
for i = 2: size(SpriN,1)
 IntSpri2(i) = trapz(newTime(1:i),SpriN(1:i));
end
ax2(2)= subplot(2,1,2)

plot(newTime,(NEU213yN-BF3yN)*5e3,newTime,IntSpri2,'r','LineWidth',2);%,newTime,IntSpri,'LineWidth',2);
grid on; 
legend('NEU213-BF3','Int(Spri)','Location','northwest'); 
xlabel('time(s)','interpreter','latex')
linkaxes(ax2,'x');
xlim([newTime(1),newTime(size(newTime,1)-1)]);

plotParam;



saveas(figure37714,strcat('plot\figure37714'), 'jpeg')
saveas(figure37714,strcat('plot\figure37714'), 'fig')
saveas(figure37714,strcat('plot\figure37714'), 'png')
saveas(figure37714,strcat('plot\figure37714'), 'epsc')
%%
% figure('Name','ok')
% [y,t,x] = lsim(tf([1],[1 0]), SpriN', newTime);
% plot(t,SpriN','k', t,y, 'r', t,IntSpri2/9.5E4, 'g-.' )
% legend('Measured', 'Simulated', 'my' )





%% suppression
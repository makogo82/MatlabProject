%% Identificazione dez
%% 6/6/1014

clc;clear all;close all;
Ts = 0.0005; 
radiusNum=32;
radiusIntPosition=1.22;

shot=38629;
shots =         [38591  39091	 39092      39093  38629  38628   38627   38404    38616 ];
shots_time =    [ 0.1    0.45     1         0.1    0.85     0.1     0.8   0.85       0.1 ];
shots_endtime = [ 1.35   0.85     1.138     1      0.95  1.195   1.22   1.12       0.6 ];

t_start = shots_time(find(shots==shot));
t_end = shots_endtime(find(shots==shot));
temp = CleanFTUdataLOCAL_v4(t_start,t_end,Ts,shot,radiusNum,radiusIntPosition);
z1 = temp.zIddData;z1Sim = temp.zIddDataSim;
plot(z1)
pause
order=[1 5 2];                      % [Ny Nu Nx]
t=0;

%%


init=[temp.dez(1); temp.dez_d_error(1)];
%init=[0; 0];
initCell = {init(1,:) init(2,:)};
initialStates = struct('Name',  {'x1' 'x2'},...
                       'Unit',  {'' ''},...
                       'Value',  initCell,...
                       'Minimum', {-Inf -Inf},...
                       'Maximum',{Inf Inf},...
                       'Fixed',{0 0});

p1 = 1E-5;
p2 = 0.8; 
p3 = 1E-8;  
p4 = 0.2;

par=[p1 p2  p3 p4];

model = idnlgrey('modelCoilH',order,par,initialStates,t,'Name', 'Model Coil H Identification');

%% parameter
setpar(model, 'Fixed', {0 0 0 0}); 
z1 = temp.zIddData;z1Sim = temp.zIddDataSim; 

%% parameter
model.Algorithm.SimulationOptions.Solver= 'ode45';
model.Algorithm.Display='On';
model.Algorithm.Tolerance=1e-15;
model.Algorithm.MaxIter=10000;
model.Algorithm.SearchMethod='Auto';
%% 
display('Merge data ');
z = merge(z1); %multiexperiment   

%identificazione del modello
display('Start Identification');
modelId = pem(z, model);
param = getpar(modelId,'Value')
%%
compare(modelId,z);

break
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


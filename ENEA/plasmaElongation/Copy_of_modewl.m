close all;
clc;
clear all;
shots =         [38609 38591  39091	39092  39094  38629  38628   38627   38404  38616 ];
shots_time =    [0.1 0.1    0.1     0.1  0.1   0.1     0.1     0.1     0.1       0.1   ];
shots_endtime = [1.1 1.35     1     1.2  1     0.9855  1.195   1.22    1.122     0.9 ];
c_dez = 2;
d_dez = 9;
c_ddez = 2;
d_ddez = 10;
for i=1:length(shots)
    shot=shots(i);
    Ts = 0.5E-3;
    samplingTime = Ts;
    t_start = shots_time(find(shots==shot));
    t_end = shots_endtime(find(shots==shot));
    data = CleanFTUdata_v6(t_start,t_end,Ts,shot,1);
    mytime = data.time;
    figure('Name',strcat('Shot', num2str(shot)));
    as(1)=subplot(4,1,1);
    plot(mytime,data.Hmis,'g',mytime,data.Hcalc,'k','LineWidth',2);
    legend('Hmis','Hcalc','Fmis','Fcalc')
    grid on; 
    as(2)=subplot(4,1,2);    
    d_error = mypseudo_derivative(mytime,data.dez,c_dez,d_dez);
    dd_error = mypseudo_derivative(mytime,d_error,c_ddez,d_ddez);
    plot(mytime,data.dez,'LineWidth',2); 
    legend('dez','Location','NorthWest')
    grid on;
    as(3)=subplot(4,1,3);
    plot(mytime,data.Fmis,mytime,data.Fcalc,'LineWidth',2);
 
    legend('d_error','Location','NorthWest')
    grid on;
    as(5)=subplot(4,1,4);
    plot(mytime,data.zs1,mytime,data.zs2,mytime,(data.zs1+data.zs2)*0.5);
    grid on; 
    title(num2str(shot));
    linkaxes(as,'x');
    
end

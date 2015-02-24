clear all; 
clc; 
close all;


shots = [ 37424 37416 37439 37420 37419 37418 37417 37416 37405 37404 37402 37401....
          37443 37429 37438 37426 37425 37424 34258 34542 34543 34546....
          34616 34640 34693 37442];

shots = [ 37416 37439 37420 37419 37418 37417 37416 37405 37404 37402 37401....
          37443 37429 37438 37426 37425 37424 37442];      
      
      
newTime = [0:0.0005:1.5]';
b=ones(1,20)/20;      
      
 for j = 1:length(shots)
    shot = shots(j);      
    dato = load(strcat('shot_', num2str(shot),'.mat'));  %dati zana   
    
    
    try
        dati = dato.Data;
        dens =  dati.esidens.y;
        denst = dati.esidens.x;
    catch exception
        dati = dato;
        try
            dens =  dati.densita.y;
            denst = dati.densita.x;
        catch exception2
            continue
        end    
    end

    NEU213y = dati.NEU213.y;
    NEU213x = dati.NEU213.x;
    BF3x = dati.BF3.x;
    BF3y = dati.BF3.y;
    IPL = dati.IPLmis.y;
    time = dati.IPLmis.x;
    Vloop = filtfilt(b,1,dati.Vloop.y);

    iplN = interp1(time,IPL,newTime);
    NEU213yFilter = filtfilt(b,1,NEU213y); 
    BF3yFilter = filtfilt(b,1,BF3y); 

    figure('Name',strcat('Shot', num2str(shot)));
    title(strcat('log_1_0 and ipl. Shot:', num2str(shot)));
    subplot(2,1,1)
    plot(newTime,iplN,'LineWidth',2); grid on;  legend('IPL')
    subplot(2,1,2)
    semilogy(NEU213x,NEU213yFilter,'k',BF3x,BF3yFilter,'g','LineWidth',2); grid on;     
    grid on; legend('NEU213','BF3','Location','northwest'); plotParam;
    xlim([newTime(1),newTime(size(newTime,1))]);
    xlabel('time(s)');

    %save(['shot_' num2str(shots(j)) '.mat'])
    %disp(['Saved shot ' num2str(shots(j)) ]);
end   



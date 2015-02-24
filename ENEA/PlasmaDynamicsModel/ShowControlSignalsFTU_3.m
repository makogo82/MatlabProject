clc
clear all
close all



shots = [35964 35965 35959 36428 36440 36569 36574 36634 36767 37681 ...
      38249	38252	38253	38254	38255	38256	38257	38258	38262	38264	38265	38266]
shots = 38249 %38616
  %t_start and t_end have to be multiple of 1ms
t_start = 0.0
t_end = 1.6
mypause = 10.0;

FindHysteresisEvidence = 0; %35965 38356
PsiDensGammaFEB = 0;
multi_runaway_1 = 0;
xFTU = 0;
IpItVloopMhd =  1;
positioning = 1;

mhd_analisis = 0;
GainScheduling = 0;
posizionamento = 0;
dipla_vloop = 0;
dipla_neutroni_misto = 0;

Ts = 0.0005; 
SMORZ1 = sqrt(1.01)
MyFontSize = 12%SMORZ1 = sqrt(1.045)

c_dipla = 2
d_dipla = 4
Soglia_dIP = 0.75E7 %(abs value)
c_dipla_slow = 4
d_dipla_slow = 26
c_dvloop = 2
d_dvloop = 4
c_dvloop_slow = 6
d_dvloop_slow = 20

c_1 = 3
d_1 = 10

c_dmhd = 10;
d_dmhd = 50;



n=1;
globalVarTokamakFTU;
%colors up to 7 shots
listaColori = ['r','b','g','k','c','y','m','r--','b--','g--','k--','c--','y--','m--','r:','b:','g:','k:','c:','y:','m:'];
MyFontSize = 14; 

for n = 1: length(shots)


%%  ELABORAZIONE DEI SEGNALI

data = CleanFTUdata_v6(t_start,t_end,Ts,shots(n),1);
mytime = data.time;
CurrentSignIp = data.CurrentSignIp;
maxIp = max(abs(data.IPLmis));


%%

Nfigure = 1;

if(xFTU==1)
    
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    clear ax
    
   
    ax(1) = subplot(4,2,1);
    plot(mytime,data.elong,'r-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Elong','Location','NorthWest','FontSize',MyFontSize)
    axis([mytime(1),mytime(end),0.95,1.23])
    
    
    mymodelDEZ = tf([0.001],[-1/10 1])*tf([1],[0.0001 1]);
    
    
    
    ax(3) = subplot(4,2,3);
    plot(mytime,data.dez,'b-',mytime, lsim(mymodelDEZ,data.Hcalc,mytime),'r--','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('dez','Location','NorthWest','FontSize',MyFontSize)
    ylim([-0.25 0.25])
    ax(5) = subplot(4,2,5);
    plot(mytime,mypseudo_derivative(mytime,data.dez,c_1,d_1),'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('\dot dez','Location','NorthWest','FontSize',MyFontSize)
    
    ax(7) = subplot(4,2,7);
    plot(mytime,data.Hmis,'b-',mytime,data.Hcalc,'r-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend(' H_m_i_s','H_c_a_l_c','Location','NorthWest','FontSize',MyFontSize)

    ax(2) = subplot(4,2,2);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','Location','NorthWest','FontSize',MyFontSize)
 
    ax(4) = subplot(4,2,4);
    plot(mytime,data.zs1,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Z1','Location','NorthWest','FontSize',MyFontSize)
    axis([mytime(1),mytime(end),0,0.335])
    
    ax(6) = subplot(4,2,6);
    plot(mytime,data.zs2,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Z2','Location','NorthWest','FontSize',MyFontSize)
    axis([mytime(1),mytime(end),-0.335,0])
    
    
    ax(8) = subplot(4,2,8);
    plot(mytime,data.dep,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend(' dep','Location','NorthWest','FontSize',MyFontSize)

    
    
    linkaxes(ax,'x');
    
    %%
    
end


if(IpItVloopMhd)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    ax(1) = subplot(3,2,1);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
    plot(mytime,data.runipprep,'r-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','I_P_p_R_E','Location','SouthEast','FontSize',MyFontSize)
    
    ax(2) = subplot(3,2,2);
    plot(mytime,mypseudo_derivative(mytime,data.IPLmis,c_1,d_1),'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('\dot I_p','Location','SouthEast','FontSize',MyFontSize)
    
    ax(3) = subplot(3,2,3);
    plot(mytime,data.Tcalc,'r-',mytime,data.Tmis,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('T_c_a_l_c','T_m_e_a_s','Location','SouthEast','FontSize',MyFontSize)
    
    ax(4) = subplot(3,2,4);
    plot(mytime,mypseudo_derivative(mytime,data.Tmis,c_1,d_1),'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('\dot T_m_i_s','Location','SouthEast','FontSize',MyFontSize)

    
    ax(5) = subplot(3,2,5);
    plot(mytime,data.vloop,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Vloop','Location','SouthEast','FontSize',MyFontSize)

    ax(6) = subplot(3,2,6);
    plot(data.mhdfstc04_time,data.mhdfstc04,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('mhdfstc04','Location','SouthEast','FontSize',MyFontSize)

    
    linkaxes(ax,'x');
    
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    ax(1) = subplot(3,2,1);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
    set(gca,'FontSize',MyFontSize)
    hold off
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','Location','SouthEast','FontSize',MyFontSize)
    
    ax(2) = subplot(3,2,2);
    plot(mytime,data.esoftx*1E-4,'b-','LineWidth',1.5);
    hold on
    plot(mytime,data.fbina_m_ris*4E14,'r-','LineWidth',1.5);
    plot(mytime,data.neutrnfc144,'k-','LineWidth',1.5);  
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Soft-X*1E-4','Hard-X*4E14','FC','Location','SouthEast','FontSize',MyFontSize)
    
    ax(3) = subplot(3,2,3);
    plot(mytime,data.rs2,'b-',mytime,data.rerif2,'r-','LineWidth',1.5);
    grid on
    ylim([0.9,1.25]);
    set(gca,'FontSize',MyFontSize)
    legend('RS2','new Rif','Location','SouthEast','FontSize',MyFontSize)
    
    ax(4) = subplot(3,2,4);
    plot(data.mhdfstc04_time,data.mhdfstc04,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('MHD','Location','SouthEast','FontSize',MyFontSize)

    
    ax(5) = subplot(3,2,5);
    
    semilogy(mytime,data.neu213,'b-',mytime,data.bf3,'r-',...
             mytime,data.densita*1E-8,'g-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('neu213','BF3','Dens*1E-8','Location','SouthEast','FontSize',MyFontSize)

    ax(6) = subplot(3,2,6);
    plot(mytime,data.vloop,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Vloop','Location','SouthEast','FontSize',MyFontSize)

    
    linkaxes(ax,'x');
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    
    ax(1) = subplot(3,2,1);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
    plot(mytime,data.runipprep,'r-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','I_P_p_R_E','Location','SouthEast','FontSize',MyFontSize)
    
    
    
    ax(2) = subplot(3,2,2);
    plot(mytime,data.Tcalc,'b:',mytime,data.Tmis,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('T_d_e_s','T_m_e_a_s','Location','SouthEast','FontSize',MyFontSize)
    
    ax(3) = subplot(3,2,3);
    plot(mytime,data.rs2,'b-',mytime,data.rerif2,'r-',...
         mytime,data.rerif1,'r:',mytime,data.rerif,'g:',...
         mytime,data.rs1,'m-', 'LineWidth',1.5);
    grid on
    ylim([0.64,1.25]);
    set(gca,'FontSize',MyFontSize)
    legend('RS2','ref2','ref1','ref','RS1','Location','SouthEast','FontSize',MyFontSize)
   
    ax(4) = subplot(3,2,4);
    plot(mytime,data.Fcalc,'b:',mytime,data.Fmis,'b-',...
         mytime,data.Vcalc,'r:',mytime,data.Vmis,'r-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('F_d_e_s','F_m_e_a_s','V_d_e_s','V_m_e_a_s','Location','SouthEast','FontSize',MyFontSize)
    
    ax(5) = subplot(3,2,5);
    plot(mytime,data.zs1,'b-',mytime,data.zurif,'r-','LineWidth',1.5);
    grid on
    ylim([0.0,0.3]);
    set(gca,'FontSize',MyFontSize)
    legend('ZS1','ref','Location','SouthEast','FontSize',MyFontSize)
    
    ax(6) = subplot(3,2,6);
    plot(mytime,data.Hcalc,'b:',mytime,data.Hmis,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('H_d_e_s','H_m_e_a_s','Location','SouthEast','FontSize',MyFontSize)
    
    
    linkaxes(ax,'x');
  
    
end
    


if(FindHysteresisEvidence)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    ax(1) = subplot(4,1,1);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','Location','SouthEast','FontSize',MyFontSize)
    %%
    ax(2) = subplot(4,1,2);
    %semilogy(mytime,data.neu213,'b-',mytime,data.bf3,'r-','LineWidth',1.5);
    plot(mytime,abs(data.neu213)./(abs(data.bf3) + 1E8),'r-','LineWidth',1.5)
    grid on
    set(gca,'FontSize',MyFontSize)
    %legend('neu213','BF3','Dens*1E-8','Location','SouthEast','FontSize',MyFontSize)
    legend('abs(neu213)/(1+abs(BF3))','Location','SouthEast','FontSize',MyFontSize)
    ylim([0 700])
    %%
    ax(3) = subplot(4,1,3);
    plot(mytime,data.vloop,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('Vloop','Location','SouthEast','FontSize',MyFontSize)
    ylim([-5 5])
    
    ax(4) = subplot(4,1,4);
    plot(data.mhdfstc04_time,abs(data.mhdfstc04*1E18),'r',mytime, data.densita,'b-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('abs(mhdfstc04)*1E18','densita','Location','SouthEast','FontSize',MyFontSize)
     ylim([0 0.5E20])
    
    linkaxes(ax,'x');
       

end


if(PsiDensGammaFEB == 1)
    
    
    psiDensityGammaLoad=load([ '/Volumes/MacintoshHD2/Uni/FTU_dati/FTU_SHOTS_DATA/Psi_Density_Gamma' num2str(shots) '.mat']);
    psiDensityGamma=psiDensityGammaLoad.Psi_Density_Gamma;
    
    %% Density
    density_x = psiDensityGamma.density.x;
    density_y = psiDensityGamma.density.data(:,1:size(density_x,2));
    sizePosMis = size(density_x,2);
    
    % interpolazione dei dati secondo la nuova posizione
    numberRadius = 32;
    Radial_Dens_position=[density_x(1):(density_x(sizePosMis)-density_x(1))/(numberRadius):density_x(sizePosMis)];
    deltaX=(Radial_Dens_position(2)-Radial_Dens_position(1))/2;
    density = max(0.1E18,interp1(density_x,density_y',Radial_Dens_position,'pchip')');
    time = psiDensityGamma.density.extrapolation(:,1);
    
    %retrieve index
    tBeginIndex = find(time >= t_start,1);
    tEndIndex = find(time >= t_end,1);
    % data scaled
    densityInterval = density(tBeginIndex:tEndIndex,:);
    timeInterval = time(tBeginIndex:tEndIndex,:);
    for iter = 1:length(Radial_Dens_position),
        LegendPos{iter} = strcat('r= ', num2str(Radial_Dens_position(iter)));
    end
    myfigure(Nfigure,0,timeInterval,densityInterval,LegendPos,'Density','time','$m^{-2}$',0)
    Nfigure = Nfigure + 1 ;
    
    %% Figure with density, and Ipl Vloop
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    ax(1) = subplot(2,1,1);
    hold on
    title(['shot #' num2str(shots(n))],'FontSize',18,'Interpreter','Latex')
    listaColori = {'r','b','g','k','c','y','m','r--','b--','g--','k--','c--','y--','m--','r:','b:','g:','k:','c:','y:','m:','r.-','b.-','g.-','k.-','c.-','y.-','m.-','r-o','b-o','g-o','k-o','c-o','y-o','m-o','r-x','b-x','g-x','k-x','c-x','y-x','m-x'};
    kk = 1;
    for iter = 1:3:length(Radial_Dens_position),
        LegendPosRed{kk} = strcat('r= ', num2str(Radial_Dens_position(iter)));
        kk = kk + 1;
        plot(timeInterval,densityInterval(:,iter),cell2mat(listaColori(kk)),'LineWidth',1.5)
    end
    hold off
    grid on
    ylabel('Interferom. $[m^{-2}]$','FontSize',16,'Interpreter','Latex')
    legend(LegendPosRed,'Location','SouthEast','FontSize',MyFontSize)
    set(gca,'FontSize',MyFontSize)
    
    ax(2) = subplot(2,1,2);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5); 
    grid on
    hold on
    plot(mytime,data.vloop*2E4,'r--','LineWidth',1.5); 
    hold off
    xlabel('time [s]','FontSize',18,'Interpreter','Latex')
    legend('I_P','V_l_o_o_p*2E4','Location','SouthEast','FontSize',MyFontSize)
    ylim([0,380E3]);
    set(gca,'FontSize',MyFontSize)
    
    linkaxes(ax,'x')
    
    %% Density, radial, time
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    t1 = t_start;
    t2 = t_end;
    j1 = find(timeInterval>=t1,1);
    j2 = find(timeInterval>=t2,1);
    hold on
    title(['shot #' num2str(shots(n))],'FontSize',18,'Interpreter','Latex')
    listaColori = {'r','b','g','k','c','y','m','r--','b--','g--','k--','c--','y--','m--','r:','b:','g:','k:','c:','y:','m:','r.-','b.-','g.-','k.-','c.-','y.-','m.-','r-o','b-o','g-o','k-o','c-o','y-o','m-o','r-x','b-x','g-x','k-x','c-x','y-x','m-x'};
    kk = 1;
    for iter = 1:1:length(Radial_Dens_position),
        kk = kk + 1;
        plot3(Radial_Dens_position(iter)*ones(1,length(timeInterval(j1:j2))),timeInterval(j1:j2),densityInterval(j1:j2,iter),cell2mat(listaColori(kk)),'LineWidth',1.5)
        hold on
    end
    hold off
    grid on
    xlabel('Major radius [m]','FontSize',MyFontSize)
    ylabel('time [s]','FontSize',MyFontSize)
    zlabel('Line density [m^{-2}]','FontSize',MyFontSize)
    legend(LegendPosRed,'Location','SouthEast','FontSize',MyFontSize)
    title(num2str(shots(n)),'FontSize',MyFontSize)
    
    
    %% Density, radial, time
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    t1 = t_start;
    t2 = t_end;
    j1 = find(timeInterval>=t1,1);
    j2 = find(timeInterval>=t2,1);
    XX = Radial_Dens_position;
    YY = timeInterval(j1:j2);
    ZZ = zeros(length(XX),length(YY));
    hold on
    title(['shot #' num2str(shots(n))],'FontSize',18,'Interpreter','Latex')
    for iter = 1:1:length(Radial_Dens_position),
        ZZ(iter,:) = densityInterval(j1:j2,iter);
    end
    surf(XX,YY,ZZ')
    grid on
    xlabel('Major radius [m]','FontSize',MyFontSize)
    ylabel('time [s]','FontSize',MyFontSize)
    zlabel('Line density [m^{-2}] ','FontSize',MyFontSize)
    title(num2str(shots(n)),'FontSize',MyFontSize)
    legend(LegendPosRed,'Location','SouthEast','FontSize',MyFontSize)
    set(gca,'FontSize',MyFontSize)
    
    
    %% FITTING THE DENSITY WITH GAUSSIAN
    DensityGaussianFit = 1;
   
    
    if(DensityGaussianFit == 1)
        %estimated complete density over all the major radius
        jG = 2; %number of gaussian terms/functions
        PausaGrafica = 0; %show the estimated curve at each time
        x1 = [rcameraInt,Radial_Dens_position];
        x_dens = [[x1(1):0.01:x1(2)] x1(3:end) ,rcameraExt];
        %tic
        GaussianFit = DensityGaussianFitFunction(shots(n),x1,x_dens,timeInterval,densityInterval,jG,PausaGrafica,1.1,1,1E18);
        %toc
    end
    
    
    %%
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    
    clear ax
    ax(1) = subplot(3,2,1);
    plot(mytime,interp1(psiDensityGamma.psi.extrapolation(:,1),psiDensityGamma.psi.extrapolation(:,2),mytime),'b--','LineWidth',1.5);
    hold on
    plot(mytime,interp1(psiDensityGamma.density.extrapolation(:,1),psiDensityGamma.density.extrapolation(:,2),mytime),'k:','LineWidth',1.5);
    plot(mytime,0.935+interp1(psiDensityGamma.gamma.extrapolation(:,1),psiDensityGamma.gamma.extrapolation(:,2),mytime),'r-','LineWidth',1.5);
    plot(mytime,0.935+0.01*interp1(psiDensityGamma.feb.extrapolation(:,1),psiDensityGamma.feb.extrapolation(:,2),mytime),'c.-','LineWidth',1.5);
    plot(mytime,data.rs2,'m-','LineWidth',1.5);
    plot(mytime,data.rs1,'m--','LineWidth',1.5);
    plot(mytime,interp1(GaussianFit.time,GaussianFit.center,mytime),'g:','LineWidth',1.5);
    
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('Magn.Center_O_D_I_N','Line Density','NC','FEB','rs2','rs1','DensGauss','Location','SouthEast','FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',MyFontSize)
    ylim([0.64,1.24]);
    ylabel('Major radius [m]','Interpreter','Latex','FontSize',12)
    grid on
    
    ax(2) = subplot(3,2,2);
    plot(mytime,interp1(psiDensityGamma.psi.extrapolation(:,1),psiDensityGamma.psi.extrapolation(:,3),mytime),'b--','LineWidth',1.5);
    hold on
    plot(mytime,1E-20*interp1(psiDensityGamma.density.extrapolation(:,1),psiDensityGamma.density.extrapolation(:,3),mytime),'k:','LineWidth',1.5);
    plot(mytime,2E-7*interp1(psiDensityGamma.gamma.extrapolation(:,1),psiDensityGamma.gamma.extrapolation(:,3),mytime),'r-','LineWidth',1.5);
    plot(mytime,1E-3*interp1(psiDensityGamma.feb.extrapolation(:,1),psiDensityGamma.feb.extrapolation(:,3),mytime),'c.-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    ylabel('Maximal values (scaled)','Interpreter','Latex','FontSize',12)
    legend('Flux_O_D_I_N','Density*1E-20 ','NC*2E-7','FEB*1E-3','Location','SouthEast','FontSize',MyFontSize)
    grid on
    
    %35965
    R(1) = find(Radial_Dens_position >= 0.88,1);
    R(2) = find(Radial_Dens_position >= 0.96303,1);
    R(3) = find(Radial_Dens_position >= 1.1453,1);
    R(4) = find(Radial_Dens_position >= 1.2087,1);
    R(5) = find(Radial_Dens_position >= 1.22,1);
    
    
    clear LegendPos;
    for iter = 1:5,
        LegendPos{iter} = strcat('r= ', num2str(Radial_Dens_position(R(iter))));
    end
    LegendPos(iter+1)={'FC*5.4E4'};
    
    ax(3) = subplot(3,2,3);
    plot(timeInterval,densityInterval(:,R(1)),'m:','LineWidth',1.5);
    hold on
    plot(timeInterval,densityInterval(:,R(2)),'g.-','LineWidth',1.5);
    plot(timeInterval,densityInterval(:,R(3)),'c--','LineWidth',1.5);
    plot(timeInterval,densityInterval(:,R(4)),'b-','LineWidth',1.5);
    plot(timeInterval,densityInterval(:,R(5)),'r-','LineWidth',1.5);
    plot(mytime,data.neutrnfc144*scaledFactorNeutrnfc144P,'k-','LineWidth',1.5);
    hold off
    grid on
    set(gca,'FontSize',MyFontSize)
    legend(LegendPos ,'Location','SouthEast','FontSize',MyFontSize)
    ylabel('$e^{-}$ number  $[m^{-2}]$','Interpreter','Latex','FontSize',12)
    
    
    %% compute range to obtain the total number of electrons
    Delta_radius = 0.005;%multiple of 1 or 5 millimiter
    new_small_radius = [R_torus:Delta_radius:R_torus+r_torus];
    y_chord = zeros(length(new_small_radius),1);%chords
    for i=1:length(new_small_radius) 
        y_chord(i)= 2*sqrt(r_torus^2-(new_small_radius(i)-R_torus).^2);
    end
    y_chord = [y_chord(end:-1:2)' y_chord'];
    
    % compute volume from the center    
    volume = zeros(length(y_chord)-1,1);
    for i=1:length(volume)
        volume(i)=2*pi*(R_torus-r_torus+(i)*Delta_radius-Delta_radius/2)*((y_chord(i)+y_chord(i+1))/2)*Delta_radius;
        %volume(i)=pi*((R_torus-r_torus+(i)*Delta_radius)^2-(R_torus-r_torus+(i-1)*Delta_radius)^2)*((y_chord(i)+y_chord(i+1))/2);
    end;
    disp('Calcolo il volume');
    sum(volume)
    errore_calcolo_volume = sum(volume) - 2*pi^2*R_torus*r_torus^2
    
    %% compute electron number
    disp('Total number of electrons');
    %Define a new quantization of the major radius
    Total_radius_quantized = [R_torus-r_torus:Delta_radius:R_torus+r_torus];
    temp1 = find(Total_radius_quantized >= Radial_Dens_position(1),1);
    temp2 = find(Total_radius_quantized <= Radial_Dens_position(end));
    reduced_radius = Total_radius_quantized(temp1:temp2(end));
    PartialElectronDensity = zeros(length(densityInterval(:,1)),length(reduced_radius));
    PartialElectronNumber = zeros(length(densityInterval(:,1)),length(reduced_radius));
    EstimatedElectronDensity = zeros(length(densityInterval(:,1)),length(Total_radius_quantized));
    EstimatedElectronIntegral = zeros(length(densityInterval(:,1)),1);
    PartialElectronIntegral = zeros(length(densityInterval(:,1)),1);
    
    for j=1:length(densityInterval(:,1))
        PartialElectronDensity(j,:) = interp1(Radial_Dens_position,densityInterval(j,:),reduced_radius)./(y_chord(temp1:temp2(end))); 
        PartialElectronNumber(j,:) = PartialElectronDensity(j,:)'.*volume(temp1:temp2(end));
        PartialElectronIntegral(j) = sum(PartialElectronNumber(j,:));
        EstimatedElectronDensity(j,:) = interp1(GaussianFit.radii,GaussianFit.mygaussian(j,:),Total_radius_quantized);
        EstimatedElectronDensity(j,2:end-1) = EstimatedElectronDensity(j,2:end-1)./(y_chord(2:end-1)); 
        EstimatedElectronIntegral(j) = sum(EstimatedElectronDensity(j,1:end-1)'.*volume);
    end;
    
    %%
    tempTime = [mytime(1):0.001:mytime(end)];
    integralFC = zeros(length(tempTime),1);
    for i = 2:length(tempTime),
        integralFC(i) = integralFC(i-1)+data.neutrnfc144(find(mytime>=tempTime(i),1));
    end
    ax(4) = subplot(3,2,4);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,interp1(timeInterval,1.5E-15*PartialElectronIntegral,mytime),'r-','LineWidth',1.5);
    plot(mytime,interp1(timeInterval,1.5E-15*EstimatedElectronIntegral,mytime),'r:','LineWidth',1.5);
    plot(mytime,data.vloop*1E4,'k-','LineWidth',1.5);
    %plot(tempTime,integralFC,'k-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('I_P','N_e*1E-14 semi','N_e*1E-14 est.','V_{loop}1E4','Location','SouthEast','FontSize',MyFontSize)
    
    
    
    %% delta plasma current, electron number, FC chamber
    DeltaT = 0.002; %It has to be a multiple of 1 ms
    tempTime = [mytime(1):DeltaT:mytime(end)];
    tempIp = interp1(mytime,data.IPLmis,tempTime);
    DeltaIp = zeros(length(tempIp),1);
    for i = 2:length(tempIp),
        DeltaIp(i) = tempIp(i)-tempIp(i-1);
    end
    %pay attention, FC is a count each 1 ms
    DeltaFC = zeros(length(DeltaIp),1);
    for i = 2:length(tempTime),
        j1 = find(mytime >= tempTime(i),1);
        DeltaFC(i) = sum(data.neutrnfc144(j1-(DeltaT/0.001)+1:j1));
    end
    DeltaNelectrons = zeros(length(DeltaIp),1);
    EstimatedDeltaNelectrons = zeros(length(DeltaIp),1);
    timeInterval(end) = mytime(end);
    for i = 2:length(tempTime),
        j1 = find(timeInterval >= tempTime(i),1);
        DeltaNelectrons(i) = PartialElectronIntegral(j1,end)-PartialElectronIntegral(j1-(DeltaT/0.0005)+1,end);
        EstimatedDeltaNelectrons(i) = EstimatedElectronIntegral(j1,end)-EstimatedElectronIntegral(j1-(DeltaT/0.0005)+1,end);
    end
    
    ax(6) = subplot(3,2,6);
    plot(tempTime,-4E10*DeltaIp,'b-','LineWidth',1.5);
    hold on
    plot(tempTime,-1E-4*DeltaNelectrons,'r-','LineWidth',1.5);
    plot(tempTime,-1E-4*EstimatedDeltaNelectrons,'r:','LineWidth',1.5);
    plot(tempTime,DeltaFC,'k-','LineWidth',1.5);
    %plot(mytime,interp1(psiDensityGamma.feb.extrapolation(:,1),psiDensityGamma.feb.extrapolation(:,4),mytime),'c.-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('-Delta I_p*4E10','-Delta N_e*1E-4','-Delta Est. N_e*1E-4','Delta FC','Location','SouthEast','FontSize',MyFontSize)
    %ylim([0.64,1.24]);
    ylabel(['Variations each ' num2str(DeltaT) 's'],'Interpreter','Latex','FontSize',12)
    grid on
    xlabel('time [s]','Interpreter','Latex','FontSize',12)
    %%
    
    
    %% Percentage of RE and their energy.....no way!
    p2 = 1 %assuming light velocity
    Bt = 6; %toroidal field
    GaussianFraction4MinorRadius = 0.7;
    
    Ne_centroid = zeros(length(tempTime),1);
    Ne_centroid_estimated = zeros(length(tempTime),1);
    Ne_minor_radius = zeros(length(tempTime),1);
    Ne_minor_radius_estimated = zeros(length(tempTime),1);
    percREonColdNe = zeros(length(tempTime),1);
    percREonColdNe_estimated = zeros(length(tempTime),1);
    theta_pitch = zeros(length(tempTime),1);
    theta_pitch_estimated = zeros(length(tempTime),1);
    
    
    for j=1:length(tempTime)
        temp1 = find(psiDensityGamma.density.extrapolation(:,1)>=tempTime(j),1);
        Ne_centroid(j) = psiDensityGamma.density.extrapolation(temp1,2);
        temp_Gaussian = psiDensityGamma.density.data(temp1,:);
        %let's search for the radius at which the Gaussian fall down
        %GaussianFraction4MinorRadius of the maximus
        [temp1 temp2] = max(temp_Gaussian);
        temp_right = find(temp_Gaussian(temp2:end)<=GaussianFraction4MinorRadius*temp1,1);
        temp_left = find(temp_Gaussian(1:temp2)>=GaussianFraction4MinorRadius*temp1,1);
        if(isempty(temp_right))
             temp_right = length(psiDensityGamma.density.x);
        end
        if(isempty(temp_left))
            temp_left = temp_right;%it can not evaluate the radius
        end
        Ne_minor_radius(j) = (-psiDensityGamma.density.x(temp_left)+psiDensityGamma.density.x(temp2-1+temp_right))/2;
        
%         plot(psiDensityGamma.density.x,temp_Gaussian,'bs')
%         title(['Time = ' num2str(tempTime(j)) 's'],'Interpreter','Latex','FontSize',MyFontSize)
%         hold on
%         plot(psiDensityGamma.density.x(temp_left),temp1,'bx')
%         plot(psiDensityGamma.density.x(temp2-1+temp_right),temp1,'bo')
%         
        
        temp_Gaussian = GaussianFit.mygaussian(find(timeInterval>=tempTime(j),1),:);
        [temp1 temp2] = max(temp_Gaussian);
        Ne_centroid_estimated(j) = x_dens(temp2);
        %let's search for the radius at which the Gaussian fall down
        %GaussianFraction4MinorRadius of the maximus
        temp_right = find(temp_Gaussian(temp2:end)<=GaussianFraction4MinorRadius*temp1,1);
        temp_left = find(temp_Gaussian(1:temp2)>=GaussianFraction4MinorRadius*temp1,1);
        if(isempty(temp_right))
             temp_right = length(x_dens)-temp2+1;
        end
        if(isempty(temp_left))
            temp_left = temp_right;%it can not evaluate the radius
        end
        Ne_minor_radius_estimated(j) = (-x_dens(temp_left)+x_dens(temp2-1+temp_right))/2;
        
        theta_pitch(j) = atan(2*pi*Ne_minor_radius(j)^2*Bt/(mu0*Ne_centroid(j)*tempIp(j)));
        percREonColdNe(j) = -(2*pi*Ne_centroid(j))./(c*p2*cos(theta_pitch(j))*e*DeltaNelectrons(j))*DeltaIp(j);
        
        
        theta_pitch_estimated(j) = atan(2*pi*Ne_minor_radius_estimated(j)^2*Bt/(mu0*Ne_centroid_estimated(j)*tempIp(j)));
        percREonColdNe_estimated(j) = -(2*pi*Ne_centroid_estimated(j))./(c*p2*cos(theta_pitch_estimated(j))*e*EstimatedDeltaNelectrons(j))*DeltaIp(j);
        
%         plot(x_dens,temp_Gaussian,'r')
%         plot(x_dens(temp_left),temp1,'rx')
%         plot(x_dens(temp2-1+temp_right),temp1,'ro')
%         ylim([1E18 4.3E19]);
%         if(tempTime(j)>=0.76)
%             pause(0.5);
%         end
%         hold off;
        
    end
    
    
    %%
    ax(5) = subplot(3,2,5);
    plot(tempTime,180/pi*theta_pitch,'b-','LineWidth',1.5);
    hold on
    plot(tempTime,180/pi*theta_pitch_estimated,':','LineWidth',1.5);
    plot(tempTime,5000*percREonColdNe,'r-','LineWidth',1.5);
    plot(tempTime,5000*min(percREonColdNe_estimated,0.1) ,'r:','LineWidth',1.5);
    plot(tempTime,-10*DeltaFC./(DeltaNelectrons.*percREonColdNe),'c.-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('\Theta pitch',' \Theta est. pitch','%N_R_E*5E3 ',' %N_R_E*5E3 est.','\Delta FC/\Delta N_R_E *10','Location','SouthEast','FontSize',MyFontSize)
    %ylim([0.64,1.24]);
    ylabel(['Evaluated each ' num2str(DeltaT) 's'],'Interpreter','Latex','FontSize',12)
    grid on
    xlabel('time [s]','Interpreter','Latex','FontSize',12)
    ylim([0 100])

    
    
    linkaxes(ax,'x');
    
    
    
    
    
    %% Density at different times
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    %Lt = [0.8 0.805  0.81 0.82 0.821 0.822 0.823];%35965
    %Lt = [0.810 0.813 0.816 0.82 0.83 0.84 ];%35965
    %Lt = [0.764 0.7645 0.765 0.7655 0.767 0.7675 ]
    %Lt = [0.8 0.81 0.82 0.83 0.832 0.833 0.834];%35964
    %Lt = [0.79 0.8 0.81 0.82 0.83 0.844 0.846];%35959
    %Lt = [0.18 0.19 0.2 0.22 0.24 0.26 ];%38519
    %Lt = [0.25 0.255 0.26 0.265 0.27 0.275 ];%38519
    
    listaColori = {'y','m','g','k','r','b','c'};
    for i=1:length(Lt),
        LtList{i} = ['t = ' num2str(Lt(i)) 's'];
        plot(Radial_Dens_position,densityInterval(find(timeInterval>=Lt(i),1),:),cell2mat(listaColori(i)),'LineWidth',1.5);
        hold on  
    end
    for i=1:length(Lt),
        plot(x_dens,GaussianFit.mygaussian(find(timeInterval>=Lt(i),1),:),[cell2mat(listaColori(i)) ':'],'LineWidth',2);
    end
    xlim([x_dens(1),x_dens(end)]);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend(LtList,'Location','NorthWest','FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',18)
    %ylim([0.64,1.24]);
    %ylabel('$N_{e^{-}}$ LOS [$m^-2$] (fitting Gaussians in dashed)','Interpreter','Latex','FontSize',16)
    %xlabel('Major radius [m]','Interpreter','Latex','FontSize',16)
    ylabel('N_{e^{-}} LOS [m^{-2}] (fitting Gaussians in dashed)','FontSize',16)
    xlabel('Major radius [m]','FontSize',16)
    grid on
    
    
   
    
    %% RE energy and q
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    clear ax
    
    
    q_mean = 2*pi*Bt*(Ne_minor_radius.^2)./(mu0*Ne_centroid.*tempIp');
    q_mean_estimated = 2*pi*Bt*(Ne_minor_radius_estimated.^2)./(mu0*Ne_centroid_estimated.*tempIp');
    centroidGamma = interp1(psiDensityGamma.gamma.extrapolation(:,1),psiDensityGamma.gamma.extrapolation(:,2),tempTime);
    q_mean_gamma = 2*pi*Bt*(Ne_minor_radius.^2)./(mu0*centroidGamma'.*tempIp');
    
    ax(1) = subplot(2,1,1);
    plot(tempTime,q_mean,'b',tempTime,q_mean_estimated,'b:',tempTime,q_mean_gamma,'r-','LineWidth',1.5)
    %xlim([x_dens(1),x_dens(end)]);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('q mean','q mean est.','q mean gamma','Location','NorthWest','FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',MyFontSize)
    %ylim([0.64,1.24]);
    ylabel('q mean ','Interpreter','Latex','FontSize',12)
    xlabel('time [s]','Interpreter','Latex','FontSize',12)
    grid on
    
end


    
if(multi_runaway_1)
    
    good = 1;
    %myshotsGood = [35634930,35962,35963,36048,36049,36060,36063,36065,36172,36173,36174,36183,36432,36477,36439,36430,35975,35676,35679];
    %myshotsBad = [35929,35959,35964,35965,36050,36059,36067,36167,36175,36182,36440,36428,36176,36179];
    
    %myshotsGood = [36569,36574,36634]
    %offset_fig = [0, 0.024, 0.004];
    %myshotsGood = [36569 35965 35959] %[1.4*1E-5, 1][1.0*1E-5, 2]
    %offset_fig = [0, -0.19, -0.175];
    %myshotsGood = [36569 35929 36428]
    %offset_fig = [0, 0.586, -0.161];
    
    
    %myshotsGood = [ 35965 ,36634]%36574]%,36634 ] %[1.4*1E-5, 1][1.0*1E-5, 2]
    %offset_fig = [ -0.19,  0.004]%0.024]%, 0.003];
    %myshotsGood = [ 35965 ,36569]%36574]%,36634 ] %[1.4*1E-5, 1][1.0*1E-5, 2]
    %offset_fig = [ -0.19,  0.001]%0.024]%, 0.003];
    %myshotsGood = [ 36428 ,36574]%36574]%,36634 ] %[1.4*1E-5, 1][1.0*1E-5, 2]
    %offset_fig = [ -0.161, 0.024]%0.024]%, 0.003];
    
    myshotsGood = [35965, 36569, 36574, 36634];
    offset_fig = [0,0.192, 0.215, 0.194];
    myshotsGood = [36569,36574,36634,35965]
    %offset_fig = [0, 0.024, 0.004,-0.19];
    %myshotsGood = [35964 35959 35929,35965]
    %offset_fig = [0, -0.175, 0.586,-0.19];
    %myshotsGood = [38513, 38519, 38515];
    %offset_fig = [0,0.01, 0.3];
    
    
    if (good==1)
        myshots = myshotsGood;
    else
        if (bad==1 )
            myshots = myshotsBad;
        end
    end

    %Carico il vettore di spari, tempo inizio e tempo di fine del fenomeno
    %runaway interessante
    temp = load('/Volumes/MacintoshHD2/Uni/FTU_dati/CQ_Psi_onset.mat');
    CQ_Psi_onset = temp.CQ_Psi_onset;
    
    clear LtList ax
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    listaColori = {'r','b','g','k','y','c','m'};    
    myxlim = [0.0 0.46];    
    for ns=1:length(myshots)
        
        jj = find(CQ_Psi_onset(:,1)==myshots(ns));
        
        LtList{ns} = ['#' num2str(myshots(ns)) ];
        
        t_i = CQ_Psi_onset(jj,2)-0.04;%PSI_magnetocentrum(1,1);
        t_f = CQ_Psi_onset(jj,3)+0.04;%PSI_magnetocentrum(end,1);
     
    
        data = CleanFTUdata_v6(t_i,t_f,Ts,myshots(ns),1);
        mytime = data.time;
        
        ax(1) = subplot(4,1,1);
        plot(mytime-offset_fig(ns),data.IPLmis,cell2mat(listaColori(ns)),'LineWidth',2);
        hold on
        grid on
        ylabel('$I_p\,$ [kA] (ref. dashed)','Interpreter','Latex','FontSize',MyFontSize);
        set(gca,'FontSize',MyFontSize)
        legend(LtList,'Location','NorthWest','FontSize',MyFontSize)
        xlim(myxlim);
        ylim([0 4E5]);
        ax(2) = subplot(4,1,2);
        plot(mytime-offset_fig(ns),data.rs2,cell2mat(listaColori(ns)),'LineWidth',2);
        hold on
        grid on
        ylabel('$R_{ext}\,$ [m] (ref. dashed)','Interpreter','Latex','FontSize',MyFontSize);
        set(gca,'FontSize',MyFontSize)
        ylim([0.93 1.25])
        xlim(myxlim);
        
        ax(3) = subplot(4,1,3);
        plot(mytime-offset_fig(ns),data.vloop,cell2mat(listaColori(ns)),'LineWidth',2);
        grid on
        hold on
        ylabel('$V_{loop}\,$ [V]','Interpreter','Latex','FontSize',MyFontSize);
        set(gca,'FontSize',MyFontSize)
        ylim([-5 20])
        xlim(myxlim);
        
        ax(4) = subplot(4,1,4);
        plot(mytime-offset_fig(ns),data.neutrnfc144,cell2mat(listaColori(ns)),'LineWidth',2);
        hold on
        grid on
        ylabel('FC (HX*2E14 dashed)','Interpreter','Latex','FontSize',MyFontSize);
        set(gca,'FontSize',MyFontSize)
        xlabel('time [s] (offset)','Interpreter','Latex','FontSize',MyFontSize);
        ylim([0 2.5E14])
        xlim(myxlim);
        
    end
    
    
    for ns=1:length(myshots)
        
        jj = find(CQ_Psi_onset(:,1)==myshots(ns));
        
        LtList{ns} = ['#' num2str(myshots(ns)) ];
        
        t_i = CQ_Psi_onset(jj,2)-0.04;%PSI_magnetocentrum(1,1);
        t_f = CQ_Psi_onset(jj,3)+0.04;%PSI_magnetocentrum(end,1);
     
    
        data = CleanFTUdata_v6(t_i,t_f,Ts,myshots(ns),1);
        mytime = data.time;
        
        ax(1) = subplot(4,1,1);
        plot(mytime-offset_fig(ns),data.runipprep,[cell2mat(listaColori(ns)) ':'],'LineWidth',2);
        xlim(myxlim);
        
        ax(2) = subplot(4,1,2);
        plot(mytime-offset_fig(ns),data.rerif2,[cell2mat(listaColori(ns)) ':'],'LineWidth',2);     
        xlim(myxlim);
        
        ax(4) = subplot(4,1,4);
        plot(mytime-offset_fig(ns),data.fbina_m_ris*2E14,[cell2mat(listaColori(ns)) ':'],'LineWidth',2);
        xlim(myxlim);
        
    end
    
    
    
    linkaxes(ax,'x');
    
    
end




if(positioning == 1)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    clear ax
    
    ax(1) = subplot(2,4,1);
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    grid on
    hold on
    if( isstruct(data.runipprep))
        plot(mytime,data.IPLmis,'k:','LineWidth',1.5);
        
    else
        plot(mytime,data.Ipprep,'b:','LineWidth',1.5);
        
    end
    
    if( isstruct(data.marteferunaway))
        plot(mytime,data.marteferunaway,'y-','LineWidth',1.5);
    end
    
    % if( ~isstruct(disruzione))
    %
    % temp = disruzione;
    % temp.y = interp1(temp.x,temp.y,mytime);
    % temp.x = mytime;
    % plot(mytime,temp.y(j1:j2)*maxIp,'c-','LineWidth',1.2);
    % end
    %
    if( ~isstruct(data.runawayplatea))
        
        temp = data.runawayplatea;
        plot(mytime,temp*maxIp,'c:','linewidth',1.5);
        
    end
    
    if( isstruct(data.marteferunaway))
        if(find(data.marteferunaway>=1,1))
            
            temp = data.IPLprep;
            j_t0 = find(temp.x>=data.marteferunaway.x(find(data.marteferunaway.y>=1,1)),1);
            for j=j_t0:length(temp.x)-1
                temp.y(j+1) = temp.y(j)/SMORZ1;
            end
            %p_fig = plot(mytime, ((temp.y(j1:j2)/max(abs(temp.y(j1:j2))) +1)),'y','LineWidth',1.2);
            p_fig = plot(mytime, ((temp.y(j1:j2)/max(abs(temp.y(j1:j2))))),'c:','LineWidth',2);
            
        end
    end
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('I_P','I_P_p','IpSMORZ1','CQ ','REplateau ','Location','SouthEast','FontSize',MyFontSize)
    hold off
    
    
    ax(5) = subplot(2,4,5);
    
    
    plot(mytime, data.Tmis,'r','LineWidth',1.2);
    
    hold on;
    grid on;
    plot(mytime, data.Tcalc,'y','LineWidth',1.2);
    
    set(gca,'FontSize',MyFontSize)
    legend('I_T_m','I_TCalc','Location','SouthEast','FontSize',MyFontSize)
    
    hold off
    
    
    ax(3) = subplot(2,4,3);
    plot(mytime,data.Vmis,'r-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Vprep,'k:','LineWidth',2);
    plot(mytime,data.Vcalc,'b:','LineWidth',1.2);
    
    %
    % if( ~isempty(AdeltaIV))
    %
    % temp = ALVcalc;
    % temp.y = interp1(temp.x,AdeltaIV.y,mytime)
    % temp.x = mytime;
    % j_t0 = find(temp.x >= t0_offset,1);
    % %offset = mean(temp.y(1:j_t0));
    % %temp.y = temp.y - offset;
    % j1 = find(temp.x >= t_i,1);
    % j2 = find(temp.x >= t_f,1);
    % plot(mytime,temp.y(j1:j2),'m-','LineWidth',2);
    %
    % end
    %
    % if( ~isempty(avbvout))
    %
    % temp = ALVcalc;
    % temp.y = interp1(temp.x,avbvout.y,mytime)
    % temp.x = mytime;
    % j_t0 = find(temp.x >= t0_offset,1);
    % offset = mean(temp.y(1:j_t0));
    % temp.y = temp.y - offset;
    % j1 = find(temp.x >= t_i,1);
    % j2 = find(temp.x >= t_f,1);
    % plot(mytime,temp.y(j1:j2),'k:','LineWidth',1.2);
    %
    % end
    
    hold off
    legend('I_V_m','I_V_p','I_V_c','AdeltaIV','HB1','Location','SouthEast','FontSize',MyFontSize)
    set(gca,'FontSize',MyFontSize)
    
    %%
    
    
    ax(2) = subplot(2,4,2);
    
    plot(mytime,data.Fmis,'r-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Fprep,'k:','LineWidth',2);
    plot(mytime,data.Fcalc,'b:','LineWidth',2);
    
    %IFff = (-data.Vmis + p_1*data.IPLmis +CurrentSigndata.IPLmis*p_0)/p_2;
    %IFff = IFff - IFff(1);
    %plot(mytime,IFff,'g-','LineWidth',2);
    
    
    if( isstruct(data.aifref))
        
        plot(mytime,data.aifref,'g:','LineWidth',2);
        
    end
    
    
    if( isstruct(data.adeltaif))
        
        plot(mytime,data.adeltaif,'m-','LineWidth',2);
        
    end
    
    legend('I_F_m','I_F_p','I_Fcalc','I_F_f_f','data.aifref','data.adeltaif','Location','SouthEast','FontSize',MyFontSize)
    
    hold off
    set(gca,'FontSize',MyFontSize)
    
    
    %%
    
    ax(6) = subplot(2,4,6);
    plot(mytime,data.dep,'r-','LineWidth',1.2);
    grid on
    legend('data.dep','Location','SouthEast','FontSize',MyFontSize)
    set(gca,'FontSize',MyFontSize)
    
    ax(7) = subplot(2,4,7);
    grid on
    hold on
    plot(mytime,data.rs1,'m-','LineWidth',1.2);
    plot(mytime,data.rs2,'c-','LineWidth',1.2);
    hold off
    axis([t_start t_end 0.6 1.25]);
    legend('data.rs1','data.rs2','Location','SouthEast','FontSize',MyFontSize)
    set(gca,'FontSize',MyFontSize)
    
    
    ax(4) = subplot(2,4,4);
    plot(mytime,data.dez,'g:','LineWidth',2);
    legend('dez','Location','SouthEast','FontSize',MyFontSize)
    grid on
    set(gca,'FontSize',MyFontSize)
    
    ax(8) = subplot(2,4,8);
    plot(mytime,data.zs1,'k-','LineWidth',2);
    grid on
    hold on
    plot(mytime,data.zs2,'k:','LineWidth',2);
    hold off
    legend('Ztop','Zdown','Location','SouthEast','FontSize',MyFontSize)
    axis([t_start t_end -0.5 0.5])
    set(gca,'FontSize',MyFontSize)
    
    
    linkaxes(ax,'x');
    
    
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    %%
    ax(1) = subplot(4,1,1);
    
    plot(mytime,data.Hprep,'k--','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Hmis,'r-','LineWidth',1.5);
    plot(mytime,data.Hcalc,'b:','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    legend('h_p_r_e_p','I_h_m','I_h_u','Location','SouthEast','FontSize',MyFontSize)
    
    ax(2) = subplot(4,1,2);
    
    plot(mytime,data.zs1,'k-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.zs2,'k--','LineWidth',1.5);
    plot(mytime,(data.zs2+data.zs1)/2.0,'r','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('data.zs1','data.zs2','Z mean','Location','SouthEast','FontSize',MyFontSize)
    
    ax(3) = subplot(4,1,3);
    
    plot(mytime,data.dez,'k-','LineWidth',1.5);
    grid on
    set(gca,'FontSize',MyFontSize)
    legend('data.dez','Location','SouthEast','FontSize',MyFontSize)
    
    
    ax(4) = subplot(4,1,4);
    
    data.Hmis_filtered = lsim(tf([1],[0.002 1]),data.Hmis,mytime);
    omega_n = 200;
    zita = 0.25;
    Ts = 0.0005;
    Kph = 0.15;%   5 *  0.22/Ts;
    Kih = 0.8; % 1.5;
    Kdh = 5E-4;
    FpidH = -(tf(Kph,1,Ts) + tf(Kih,[1 -1],Ts));
    Fdata.dez = -0.21*tf([1 0],[0.1 1])*tf([1],[1/omega_n^2 2*zita/omega_n 1]);
    [A,B,C,D] = ssdata(c2d(Fdata.dez,Ts));
    
    data.dez_sym = lsim(Fdata.dez,data.Hcalc,mytime);
    u_sym = lsim(FpidH,data.dez,mytime);
    integralState  = (data.Hcalc(1)-Kph*data.dez(1)/(-1.5E-4) )/Kih;
    previousError = 0;
    samplingTime = Ts;
    derivativeState = 0;
    tau2 = 2.0E-3;
    
    data.dez_sym(1) = data.dez(1);
    O = obsv(A,C);
    X = pinv(O)*data.dez(1:3)'/1000;
    closed_loop = 0;
    
    for k = 2:length(mytime)
        
        %if(closed_loop)
        %    error =  data.dez_sym(k-1)/(-1.5E-4);
        %else
        error = data.dez(k-1)/(-1.5E-4);
        %end
        
        proportionalBranch = Kph * error;
        
        integralState = integralState +(error + previousError) * samplingTime *0.5;
        integralBranch = Kih * integralState;
        
        derivativeState = ((error - previousError) + tau2*derivativeState) / (tau2 + samplingTime);
        derivativeBranch = Kdh * derivativeState;
        
        previousError=error;
        
        u_sym(k) = proportionalBranch + integralBranch + derivativeBranch;
        
        %if(closed_loop)
        %u = u_sym(k);
        %else
        u = data.Hcalc(k);
        %end
        
        X = A*X + B*u;
        data.dez_sym(k)= C*X+D*u;
        
    end
    
    plot(mytime,data.Hmis_filtered,'r-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Hcalc,'b:','LineWidth',1.5);
    plot(mytime,3000*data.dez,'k-','LineWidth',1.5);
    plot(mytime,3000*data.dez_sym,'c-','LineWidth',1.5);
    plot(mytime,3000*(data.zs2+data.zs1)/2.0,'g','LineWidth',1.5);
    plot(mytime,u_sym,'m-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('I_h_m','I_h_u','data.dez','mydata.dez','Z mean','my_u','Location','SouthEast','FontSize',MyFontSize)
    
    
    
    
    
    
    plot(mytime,data.Hmis_filtered,'r-','LineWidth',1.5);
    grid on
    hold on
    plot(mytime,data.Hcalc,'b:','LineWidth',1.5);
    plot(mytime,3000*data.dez,'k-','LineWidth',1.5);
    plot(mytime,data.dez_sym,'c-','LineWidth',1.5);
    plot(mytime,3000*(data.zs2+data.zs1)/2.0,'g','LineWidth',1.5);
    plot(mytime,u_sym,'m-','LineWidth',1.5);
    hold off
    set(gca,'FontSize',MyFontSize)
    legend('I_h_m','I_h_u','data.dez','mydata.dez','Z mean','my_u','Location','SouthEast','FontSize',MyFontSize)
    
    linkaxes(ax,'x');
    clear ax;
    %%
end





if(mhd_analisis == 1)
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    ax(1) = subplot(3,1,1);
    
    title(['shot #' num2str(shots(n))],'FontSize',16)
    
    ddata.IPLmis = mypseudo_derivative(mytime,data.IPLmis,c_dipla,d_dipla);
    plot(mytime,ddata.IPLmis,'g-','LineWidth',1.5);
    hold on
    ddata.IPLmis = mypseudo_derivative(mytime,data.IPLmis,c_dipla_slow,d_dipla_slow);
    plot(mytime,ddata.IPLmis,'k-','LineWidth',1.3);
    hold off
    grid on
    
    legend('d I_P/dt ', 'd I_P/dt slow ','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    %%
    ax(2) = subplot(3,1,2); 
    
    plot(data.mhdfstc04_time,data.mhdfstc04,'b-','LineWidth',1.5);
    hold on
    s1 = 0.01;
    s2 = 0.7;
    tau = 0.9995;
    n1 = 4;
    n2 = 400;
    x0 = 0;
    my_amp_mhd = zeros(1,length(data.mhdfstc04_time)); 
    for i=1:length(data.mhdfstc04_time)
       my_amp_mhd(i) = my_hysteresis_online_2r(abs(data.mhdfstc04(i)),s1,s2,n1,n2,tau,x0,0);
    end
    plot(data.mhdfstc04_time,my_amp_mhd,'r-','LineWidth',1.5);
    hold off
    grid on
    legend('mhdftsch04','my_amp','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
 
    
    
    c_dmhd_slow = 600;
    d_dmhd_slow = 2500;

    ax(3) = subplot(3,1,3);
    hold on
    %d_my_mhd_amp = mypseudo_derivative(mhdfstc04_time,my_amp_mhd,c_dmhd,d_dmhd);
    %plot(mhdfstc04_time,d_my_mhd_amp,'m-','LineWidth',1.5);
    d_my_mhd_amp = mypseudo_derivative(data.mhdfstc04_time,my_amp_mhd,c_dmhd_slow,d_dmhd_slow);
    plot(data.mhdfstc04_time,d_my_mhd_amp,'k-','LineWidth',1.5);
    hold off
    grid on
    legend('d mhd/dt slow','d mhd/dt','Location','SouthEast')
    %%
    
%     ax(5) = subplot(5,1,5);
%     
%     temp = BF3;
%     temp.y = mypseudo_derivative(temp.x,abs(temp.y),c_bf3,d_bf3);
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     if(isempty(j2))
%         j2 = length(temp.x);
%     end
%     plot(mytime,temp.y(j1:j2)/max(temp.y(j1:j2)),'b','LineWidth',1.5);
%     grid on
%     
%     legend('dBF3avrat/dt','Location','SouthEast')
%     
%     
    set(gca,'FontSize',MyFontSize)
    linkaxes(ax,'x');clear ax
    
    %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_d'])
    %%
end



%%


if(posizionamento==1)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    hold on
    temp = data.dep;
    j_t0 = find(temp.x >= t0_offset,1);
    offset = mean(temp.y(1:j_t0));
    temp.y = temp.y - offset;
    j1 = find(temp.x >= t_i,1);
    j2 = find(temp.x >= t_f,1);
    plot(mytime,temp.y(j1:j2)*1E4,'r-','LineWidth',1.2);
    
    
    temp = ALFcalc;
    j_t0 = find(temp.x >= t0_offset,1);
    offset = mean(temp.y(1:j_t0));
    temp.y = temp.y - offset;
    j1 = find(temp.x >= t_i,1);
    j2 = find(temp.x >= t_f,1);
    plot(mytime,temp.y(j1:j2),'b:','LineWidth',2);
    
    %plot(mytime, ALFcalc.y(j1) + (0.015*(data.IPLmisLmis.y(j1)-data.IPLmisLmis.y(j1:j2))),'k:','LineWidth',2);
    temp = ALVmis;
    j_t0 = find(temp.x >= t0_offset,1);
    offset = mean(temp.y(1:j_t0));
    temp.y = (temp.y - offset)*CurrentSign;
    %plot(mytime(j1:j2), -0.024*data.IPLmisLmis.y(j1:j2)-1.3*interp1(temp.x,temp.y,mytime(j1:j2)),'k:','LineWidth',2);
    plot(mytime(j1:j2), -0.045*data.IPLmisLmis.y(j1:j2)-4.0*interp1(temp.x,temp.y,mytime(j1:j2)),'k:','LineWidth',2);
    plot(mytime(j1:j2), -a1*data.IPLmisLmis.y(j1:j2)-272/63*interp1(temp.x,temp.y,mytime(j1:j2)),'y:','LineWidth',2);
    
    plot(mytime(j1:j2), interp1(ALdata.Fprep.x,ALdata.Fprep.y,mytime(j1:j2)),'g:','LineWidth',2);
    
    grid on
    
    
    legend('data.dep','F_calc','I_F new prep','Old I_f prep','I_F prep','Location','SouthEast')
    
    hold off
    set(gca,'FontSize',MyFontSize)
    
    
    %%
end




%%

if(dipla_vloop==1)
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    ax(1) = subplot(3,1,1);
    
    title(['shot #' num2str(shots(n))],'FontSize',16)
    
    ddata.IPLmis = mypseudo_derivative(mytime,data.IPLmis,c_dipla,d_dipla);
    plot(mytime,ddata.IPLmis,'g-','LineWidth',1.5);
    hold on
    ddata.IPLmis = mypseudo_derivative(mytime,data.IPLmis,c_dipla_slow,d_dipla_slow);
    plot(mytime,ddata.IPLmis,'k-','LineWidth',1.3);
    hold off
    grid on
    
    legend('d I_P/dt ', 'd I_P/dt slow ','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    
    ax(2) = subplot(4,1,2)
    
    plot(mytime,data.vloop,'b-','LineWidth',1.5);
    hold on
    
    
    temp = mypseudo_derivative(mytime,data.vloop,c_dvloop,d_dvloop);
    plot(mytime,max(data.vloop)*temp/max(temp),'g:','LineWidth',1.5);
    
    temp = mypseudo_derivative(mytime,data.vloop,c_dvloop_slow,d_dvloop_slow);
    plot(mytime,max(data.vloop)*temp/max(temp),'k-','LineWidth',1.2);
    
    hold off
    
    legend('Vloop','dVloop/dt','dVloop/dt slow','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    grid on
    
    ax(3) = subplot(4,1,3)
    plot(mytime,data.densita,'r-','LineWidth',1.2); 
    legend('Densita','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    grid on
    
    
    ax(4) = subplot(4,1,4)
    plot(mytime,data.avalvecmd,'c-','LineWidth',1.2);
    
    
    legend('Deuterio','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    grid on
    
    linkaxes(ax,'x');clear ax
    
    %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_b'])
end






%%
if(dipla_neutroni_misto == 1)
    
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1 ;
    
    ax(1) = subplot(5,1,1);
    
    plot(mytime,data.IPLmis,'b-','LineWidth',1.5);
    set(gca,'FontSize',MyFontSize)
    title(['shot #' num2str(shots(n))],'FontSize',16)
    
    grid on
    
    ax(2) = subplot(5,1,2);
    
    plot(mytime,data.fbina_m_ris,'k-','LineWidth',1.2);
    grid on
    legend('Fbina m ris','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    
    ax(3) = subplot(5,1,3);
    
    plot(mytime,10*data.neu213,'b-','LineWidth',1.2);
    grid on
    hold on
    plot(mytime,50*data.bf3,'r-','LineWidth',1.2);
    plot(mytime,data.neutrnfc144,'g','LineWidth',1.2);
    hold off
    
    legend('NEU213*10','BF3*50','fc144','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    
    ax(4) = subplot(5,1,4);
    
    plot(data.mhdfstc04_time,data.mhdfstc04,'y-','LineWidth',1.2);
    grid on
    legend('mhdftsch04.ch','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    
    ax(5) = subplot(5,1,5);
    
    plot(mytime,data.temperatura,'r-','LineWidth',1.2);
    grid on
    legend('Temperatura','Location','SouthEast')
    set(gca,'FontSize',MyFontSize)
    
    
    linkaxes(ax,'x');clear ax
    %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_c'])
    
    %%
end


 
% if(vertical_field)
%     %%
%     
%     figure(Nfigure)
%     Nfigure = Nfigure + 1 ;
%     
%     ax(1) = subplot(4,1,1)
%     
%     temp = data.IPLmisLmis;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'b-','LineWidth',1.5);
%     grid on
%     legend('I_p','Location','SouthEast')
%     set(gca,'FontSize',MyFontSize)
%     title(['shot #' num2str(shots(n))],'FontSize',16)
%     
%     
%     ax(2) = subplot(4,1,2)
%     
%     temp1 = ALVmis;
%     j_t0 = find(temp1.x >= t0_offset,1);
%     offset = mean(temp1.y(1:j_t0));
%     temp1.y = (temp1.y - offset)*CurrentSigndata.IPLmis;
%     temp1.y = interp1(temp1.x,temp1.y,mytime);
%     
%     temp2 = ALFmis;
%     j_t0 = find(temp2.x >= t0_offset,1);
%     offset = mean(temp2.y(1:j_t0));
%     temp2.y = temp2.y - offset;
%     temp2.y = interp1(temp2.x,temp2.y,mytime);
%     
%     plot(mytime,alpha*temp1.y(j1:j2) + temp2.y(j1:j2) + a1*data.IPLmisLmis.y(j1:j2) + a0,'r-','LineWidth',1.5);
%     %  hold on
%     %  plot(mytime,a_V*temp1.y(j1:j2) + temp2.y(j1:j2) + a_data.IPLmis*data.IPLmisLmis.y(j1:j2) + a0,'k-','LineWidth',1.5);
%     %  hold off
%     % grid on
%     grid on
%     set(gca,'FontSize',MyFontSize)
%     legend('Horizontal Balance Equation','No Alpha Vert. Bal.','Location','SouthEast')
%     
%     
%     % ax(3) = subplot(4,1,3)
%     %
%     % temp2 = ALFmis;
%     % j_t0 = find(temp2.x >= t0_offset,1);
%     % offset = mean(temp2.y(1:j_t0));
%     % temp2.y = temp2.y - offset;
%     % temp2.y = interp1(temp2.x,temp2.y,mytime);
%     %
%     % plot(mytime,temp1.y(j1:j2),mytime,temp2.y(j1:j2),'LineWidth',1.5);
%     % grid on
%     % set(gca,'FontSize',MyFontSize)
%     % legend('I_V','I_F','Location','SouthEast')
%     
%     
%     
%     ax(3) = subplot(4,1,3)
%     
%     temp = densita
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'g-','LineWidth',1.2);
%     
%     legend('Density','Location','SouthEast')
%     set(gca,'FontSize',MyFontSize)
%     grid on
%     
%     
%     
%     
%     ax(4) = subplot(4,1,4)
%     
%     for jj = j1:j2
%         if(densita.y(jj)>=21)
%             new_a1(1+jj-j1) = a_data.IPLmis + ((a1-a_data.IPLmis)/(45.5-29.8)) * (densita.y(jj)-29.8);
%         else
%             new_a1(1+jj-j1) = a1;
%         end
%     end
%     plot(mytime,alpha*temp1.y(j1:j2) + temp2.y(j1:j2) + a1*data.IPLmisLmis.y(j1:j2) + a0,'r-','LineWidth',1.5);
%     hold on
%     plot(mytime,alpha*temp1.y(j1:j2) + temp2.y(j1:j2) + new_a1.*data.IPLmisLmis.y(j1:j2) + a0,'k-','LineWidth',1.5);
%     hold off
%     grid on
%     set(gca,'FontSize',MyFontSize)
%     legend('Horizontal Balance Equation','HBE with density corr.','Location','SouthEast')
%     
%     
%     
%     
%     linkaxes(ax,'x')
%     clear ax
%     
%     
%     %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_e'])
%     
% end
% 
% 
% if(derivata_data.IPLmis==1)
%     %%
%     
%     figure(Nfigure)
%     Nfigure = Nfigure + 1 ;
%     
%     
%     ax(1) = subplot(4,1,1)
%     
%     temp = data.IPLmisLmis;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'b-','LineWidth',1.3);
%     title(['shot #' num2str(shots(n))],'FontSize',16)
%     grid on
%     
%     hold on
%     temp = data.marteferunaway;
%     temp.y = interp1(temp.x,temp.y,mytime)
%     temp.x = mytime;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,(temp.y(j1:j2)/max(abs(temp.y(j1:j2))))*3E5,'y-','LineWidth',1.3);
%     
%     
%     temp = fbcordata.IPLmis1;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'m-','LineWidth',1.3);
%     
%     temp = runawayplatea;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)*3E5,'r-','LineWidth',1.3);
%     
%     
%     temp = disruzione;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)*3E5,'g-','LineWidth',1.3);
%     
%     
%     set(gca,'FontSize',MyFontSize)
%     legend('I_p','Runaway','data.IPLmis RE','RE plateau','disruzione','Location','SouthEast')
%     hold off
%     %%
%     ax(2) = subplot(4,1,2)
%     
%     temp = data.IPLmisLmis;
%     temp.y = mypseudo_derivative(mytime,data.IPLmisLmis.y,c_dipla,d_dipla);
%     j_t0 = find(temp.x >= t0_offset,1);
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'g-','LineWidth',1.5);
%     grid on
%     hold on
%     plot(mytime,-Soglia_ddata.IPLmis*ones(1,length(mytime)),'r-',...
%         mytime,Soglia_ddata.IPLmis*ones(1,length(mytime)),'r','LineWidth',1.5);
%     hold off
%     set(gca,'FontSize',MyFontSize)
%     legend('dI_p/dt','Soglia','Location','SouthEast')
%     %%
%     ax(3) = subplot(4,1,3)
%     
%     temp = data.IPLmisLmis;
%     temp.y = mypseudo_derivative(mytime,fbcordata.IPLmis1.y,c_dipla_slow,d_dipla_slow) - mypseudo_derivative(mytime,data.IPLmisLmis.y,c_dipla_slow,d_dipla_slow);
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'k-','LineWidth',1.5);
%     grid on
%     set(gca,'FontSize',MyFontSize)
%     legend('dI_p/dt slow','Location','SouthEast')
%     
%     
%     
%     
%     ax(4) = subplot(4,1,4)
%     mytemp = ones(1,j2-j1+1);
%     a11 = 0.8
%     for jj=1:length(mytemp)-1
%         dzSignal = 0.0;
%         if(temp.y(jj+j1-1) >= 0.5E6)
%             dzSignal = abs(temp.y(jj+j1-1))*0.5E-6;
%         end
%         mytemp(jj+1) = a11*mytemp(jj) +  (1-a11) * ( 1./(1+0.3*dzSignal^2) );
%     end
%     plot(mytime,(temp.y(j1:j2)./(1E7)),'k-','LineWidth',1.5);
%     hold on
%     plot(mytime,mytemp,'r-','LineWidth',1.5);
%     hold off
%     grid on
%     set(gca,'FontSize',MyFontSize)
%     legend('dVB1 sat_gain','Location','SouthEast')
%     
%     
%     
%     linkaxes(ax,'x')
%     clear ax
%     
%     %%
% end
% 
% 
% 
% 
% if(interazione_parete_xduri == 1)
%     
%     %%
%     figure(Nfigure)
%     Nfigure = Nfigure + 1 ;
%     
%     ax(1) = subplot(3,1,1)
%     
%     temp = data.IPLmisLmis;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)/max(temp.y(j1:j2)),'b-','LineWidth',1.5);
%     
%     grid on
%     
%     hold on
%     temp = data.marteferunaway;
%     temp.y = interp1(temp.x,temp.y,mytime)
%     temp.x = mytime;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)/max(abs(temp.y(j1:j2))),'y-','LineWidth',1.5);
%     
%     
%     temp = data.IPLmisLmis;
%     temp.y = mypseudo_derivative(mytime,data.IPLmisLmis.y,c_dipla,d_dipla);
%     j_t0 = find(temp.x >= t0_offset,1);
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)/max(abs(temp.y(j1:j2))),'g:','LineWidth',1.5);
%     
%     temp = data.IPLmisLmis;
%     temp.y = mypseudo_derivative(mytime,data.IPLmisLmis.y,c_dipla_slow,d_dipla_slow);
%     j_t0 = find(temp.x >= t0_offset,1);
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)/max(abs(temp.y(j1:j2))),'k-','LineWidth',1.3);
%     
%     
%     set(gca,'FontSize',MyFontSize)
%     title(['shot #' num2str(shots(n))],'FontSize',16)
%     legend('I_P','Prot12RE','ddata.IPLmis/dt','ddata.IPLmis/dt slow','Location','SouthEast')
%     
%     hold off
%     
%     
%     ax(2) = subplot(3,1,2)
%     
%     temp = NEU213;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(temp.x(j1:end),temp.y(j1:end)./max(temp.y(j1:end)),'b:','LineWidth',1.3);
%     grid on
%     hold on
%     temp = FBINA_M_RIS;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'k-','LineWidth',1.5);
%     
%     hold off
%     
%     legend('NEU213',temp.ch,'Location','SouthEast')
%     set(gca,'FontSize',MyFontSize)
%     
%     %%
%     
%     ax(3) = subplot(3,1,3)
%     
%     
%     
%     
%     % temp = data.dep;
%     % j_t0 = find(temp.x >= t0_offset,1);
%     % offset = mean(temp.y(1:j_t0));
%     % temp.y = temp.y - offset;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,temp.y(j1:j2)/max(abs(temp.y(j1:j2))),'r-','LineWidth',1.2);
%     %
%     
%     temp = data.data.rs1;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     %plot(mytime,30*(temp.y(j1:j2)-0.7),'k-','LineWidth',1.2);
%     plot(mytime,temp.y(j1:j2),'k-','LineWidth',1.2);
%     grid on
%     
%     hold on
%     
%     temp = data.data.rs2;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     %plot(mytime,30*(temp.y(j1:j2)-1.2),'r-','LineWidth',1.2);
%     plot(mytime,temp.y(j1:j2),'r-','LineWidth',1.2);
%     
%     
%     temp = rerif1;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     %plot(mytime,30*(temp.y(j1:j2)-1.2),'r-','LineWidth',1.2);
%     plot(mytime,temp.y(j1:j2),'b-','LineWidth',1.2);
%     
%     
%     %hold off
%     
%     % temp = data.zs1;
%     % j_t0 = find(temp.x >= t0_offset,1);
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,temp.y(j1:j2)/2,'k-','LineWidth',1.2);
%     %
%     %
%     % temp = data.zs2;
%     % j_t0 = fi4d(tempx >= t0_offset,1);
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,temp.y(j1:j2),'k:','LineWidth',2);
%     
%     %legend('data.dep','data.data.rs1','data.data.rs2','Ztop','Zdown','Location','SouthEast')
%     
%     legend('Ri','Re','Re_p_r_e_p','Location','SouthEast')
%     
%     grid on
%     axis([t_start t_end 0.6 1.3])
%     
%     set(gca,'FontSize',MyFontSize)
%     
%     %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_a'])
%     linkaxes(ax,'x');clear ax
%     
%     
%     
%     
%     figure(Nfigure)
%     
%     
%     ax(1) = subplot(4,1,1)
%     
%     temp = data.IPLmisLmis;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2)/max(temp.y(j1:j2)),'b-','LineWidth',1.5);
%     
%     grid on
%     
%     set(gca,'FontSize',MyFontSize)
%     title(['shot #' num2str(shots(n))],'FontSize',16)
%     legend('I_P','Location','SouthEast')
%     
%     hold off
%     
%     
%     ax(2) = subplot(4,1,2)
%     
%     temp = NEU213;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(temp.x(j1:end),temp.y(j1:end)./max(temp.y(j1:end)),'b:','LineWidth',1.3);
%     grid on
%     hold on
%     temp = FBINA_M_RIS;
%     j_t0 = find(temp.x >= t0_offset,1);
%     offset = mean(temp.y(1:j_t0));
%     temp.y = temp.y - offset;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,temp.y(j1:j2),'b-','LineWidth',1.5);
%     
%     hold off
%     legend('neu213','Fbina','Location','SouthEast')
%     
%     ax(3) = subplot(4,1,3)
%     
%     temp = data.data.rs1;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,(temp.y(j1:j2)+data.data.rs2.y(j1:j2))/2,'b-','LineWidth',1.2);
%     grid on
%     legend('RE bari','Location','SouthEast')
%     
%     ax(4) = subplot(4,1,4)
%     
%     
%     
%     temp = data.zs1;
%     j1 = find(temp.x >= t_start,1);
%     j2 = find(temp.x >= t_end,1);
%     plot(mytime,(temp.y(j1:j2)+data.zs2.y(j1:j2))/2,'b-','LineWidth',1.2);
%     
%     grid on
%     legend('Z bari','Location','SouthEast')
%     
%     %
%     % load('/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/shot_36183.mat');
%     %
%     %
%     % ax(1) = subplot(4,1,1)
%     % hold on
%     % temp = data.IPLmisLmis;
%     % j_t0 = find(temp.x >= t0_offset,1);
%     % offset = mean(temp.y(1:j_t0));
%     % temp.y = temp.y - offset;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,temp.y(j1:j2)/max(temp.y(j1:j2)),'r-','LineWidth',1.5);
%     %
%     % grid on
%     %
%     % set(gca,'FontSize',MyFontSize)
%     % title('I_P','FontSize',16)
%     % legend('82','83','Location','SouthEast')
%     %
%     % hold off
%     %
%     %
%     % ax(2) = subplot(4,1,2)
%     % hold on
%     %
%     % temp = NEU213;
%     % j_t0 = find(temp.x >= t0_offset,1);
%     % offset = mean(temp.y(1:j_t0));
%     % temp.y = temp.y - offset;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(temp.x(j1:end),temp.y(j1:end)./max(temp.y(j1:end)),'r:','LineWidth',1.3);
%     % grid on
%     % temp = FBINA_M_RIS;
%     % j_t0 = find(temp.x >= t0_offset,1);
%     % offset = mean(temp.y(1:j_t0));
%     % temp.y = temp.y - offset;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,temp.y(j1:j2),'r-','LineWidth',1.5);
%     %
%     % hold off
%     % title('Neutroni e X duri','FontSize',16)
%     %
%     % legend('NEU_8_2','X_8_2','NEU_8_3','X_8_3','Location','SouthEast')
%     % set(gca,'FontSize',MyFontSize)
%     %
%     %
%     % ax(3) = subplot(4,1,3)
%     %
%     % hold on
%     %
%     % temp = data.data.rs1;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,(temp.y(j1:j2)+data.data.rs2.y(j1:j2))/2,'r-','LineWidth',1.2);
%     % hold off
%     %
%     % set(gca,'FontSize',MyFontSize)
%     % title('Baricentro Orizzontale','FontSize',16)
%     % legend('RS_8_2','RS_8_3','Location','SouthEast')
%     %
%     %
%     % ax(4) = subplot(4,1,4)
%     % hold on
%     %
%     %
%     % temp = data.zs1;
%     % j1 = find(temp.x >= t_i,1);
%     % j2 = find(temp.x >= t_f,1);
%     % plot(mytime,(temp.y(j1:j2)+data.zs2.y(j1:j2))/2,'r-','LineWidth',1.2);
%     %
%     % grid on
%     % set(gca,'FontSize',MyFontSize)
%     % title('Baricentro Verticale','FontSize',16)
%     % legend('RZ_8_2','RZ_8_3','Location','SouthEast')
%     %
%     % %print('-dpng','-r400', ['Bshot_' num2str(shots(n)) '_a'])
%     linkaxes(ax,'x');
%     
%     %
%     % %%
% end
% 
% %%
% if(interazione_parete_xduri2 == 1)
%     
%     
%     colori = ['k','b','r','g','c'];
%     
%     for jj = 1:5
%         
%         
%         load(['/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/shot_' num2str(shots(end-jj+1)) '.mat']);
%         clear i j k
%         
%         Nfigure = Nfigure +1;
%         figure(Nfigure)
%         hold on
%         
%         ax(1) = subplot(4,1,1)
%         hold on
%         temp = data.IPLmisLmis;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(mytime,temp.y(j1:j2),colori(jj),'LineWidth',1.5);
%         grid on
%         
%         set(gca,'FontSize',MyFontSize)
%         title(['shot #' num2str(shots(n))],'FontSize',16)
%         legend('I_P','Location','SouthEast')
%         
%         
%         ax(2) = subplot(4,1,2)
%         hold on
%         temp = Vloop;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         %temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(mytime,lsim(tf(1,[1 0]),temp.y(j1:j2),mytime),colori(jj),'LineWidth',1.5);
%         legend('V loop Integral','Location','SouthEast')
%         %plot(mytime,lsim(tf(1,[1 0]),temp.y(j1:j2).*interp1(mytime,data.IPLmisLmis.y,mytime),mytime),colori(jj),'LineWidth',1.5);
%         %legend('V loop*I Integral','Location','SouthEast')
%         grid on
%         
%         set(gca,'FontSize',MyFontSize)
%         title(['shot #' num2str(shots(n))],'FontSize',16)
%         
%         
%         hold off
%         
%         ax(3) = subplot(4,1,3)
%         hold on
%         temp = data.data.rs2;
%         %temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(mytime,(data.data.rs2.y(j1:j2)+data.data.rs1.y(j1:j2))/2.0,colori(jj),'LineWidth',1.5);
%         plot(mytime,(data.zs2.y(j1:j2)+data.zs1.y(j1:j2))/2.0 + 0.9,[colori(jj) ':'],'LineWidth',1.5);
%         grid on
%         
%         set(gca,'FontSize',MyFontSize)
%         title(['shot #' num2str(shots(n))],'FontSize',16)
%         legend('RE baricentro','Z baricentro +0.9m','Location','SouthEast')
%         
%         hold off
%         
%         ax(4) = subplot(4,1,4)
%         hold on
%         temp = NEU213;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(temp.x(j1:end),temp.y(j1:end)./max(temp.y(j1:end)),[colori(jj) ':'],'LineWidth',1.3);
%         grid on
%         hold on
%         temp = FBINA_M_RIS;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(mytime,temp.y(j1:j2),colori(jj),'LineWidth',1.5);
%         
%         hold off
%         legend('neu213','Fbina','Location','SouthEast')
%         
%         hold off
%         
%         linkaxes(ax,'x');
%         
%         
%     end
%     
%     Nfigure = Nfigure +1;
%     figure(Nfigure)
%     for jj = 1:5
%         
%         load(['/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/shot_' num2str(shots(end-jj+1)) '.mat']);
%         clear i j k
%         
%         hold on
%         temp = NEU213;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(temp.x(j1:end),temp.y(j1:end),[colori(jj) ':'],'LineWidth',1.3);
%     end
%     
%     
%     
%     
% end
% 
% 
% 
% 
% %%
% if(baricentri_PSI == 1)
%     
%     
%     load('/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/CQ_onset.mat');
%     Nfigure = Nfigure +1;
%     
%     if(find(CQ_onset(:,1)==shots(n)))
%     %for jj = 1:length(CQ_onset(:,1))
%         jj = find(CQ_onset(:,1)==shots(n));
%         figure(Nfigure)
%         
%         load(['/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/PSI_magnetocentrum_' num2str(CQ_onset(jj,1)) '.mat']);
%         load(['/Volumes/MacintoshDati/Uni/FTU_dati/FTU_SHOTS_DATA/shot_' num2str(CQ_onset(jj,1)) '.mat']);
%         clear i j k
%         
%         t_start = PSI_magnetocentrum(1,1);
%         t_end = PSI_magnetocentrum(end,1);
%         n_psi = length(PSI_magnetocentrum(:,1));
%         
%         newtimepsi = [t_start:0.002:t_end];
%         
%         ax(1) = subplot(5,1,1);
%         hold on
%         temp = data.IPLmisLmis;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         plot(newtimepsi,interp1(temp.x,temp.y,newtimepsi),'b','LineWidth',1.5);
%         grid on
%         
%         set(gca,'FontSize',MyFontSize)
%         title(['shot #' num2str(shots(n))],'FontSize',16)
%         legend('I_P','Location','SouthEast')
%         
%         
%         ax(2) = subplot(5,1,2);
%         hold on
%         temp = Vloop;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         %temp.y = temp.y - offset;
%         plot(newtimepsi,interp1(temp.x,temp.y,newtimepsi),'b','LineWidth',1.5);
%         legend('V loop Integral','Location','SouthEast')
%         %plot(mytime,lsim(tf(1,[1 0]),temp.y(j1:j2).*interp1(mytime,data.IPLmisLmis.y,mytime),mytime),colori(jj),'LineWidth',1.5);
%         %legend('V loop*I Integral','Location','SouthEast')
%         grid on
%         
%         set(gca,'FontSize',MyFontSize)
%         hold off
%         
%         
%         
%         ax(3) = subplot(5,1,3);
%         hold on
%         temp = data.data.rs2;
%         %temp.y = temp.y - offset;
%         geo_bar = (interp1(data.data.rs2.x,data.data.rs2.y,newtimepsi)+interp1(data.data.rs1.x,data.data.rs1.y,newtimepsi))/2.0;
%         plot(newtimepsi,geo_bar,'b','LineWidth',1.5);
%         %plot(newtimepsi,(interp1(data.zs2.x,data.zs2.y,newtimepsi)+interp1(data.zs1.x,data.zs1.y,newtimepsi))/2.0,'b:','LineWidth',1.5);
%         plot(newtimepsi,PSI_magnetocentrum(:,2),'r','LineWidth',1.5);
%         
%         %plot(mytime,(data.data.data.rs2.y(j1:j2)+data.data.rs1.y(j1:j2))/2.0,'r:','LineWidth',1.5);
%         %plot(mytime,(data.zs2.y(j1:j2)+data.zs1.y(j1:j2))/2.0 + 0.9,'r:','LineWidth',1.5);
%         grid on
%         set(gca,'FontSize',MyFontSize)
%         %legend('RE baricentro','Z baricentro +0.9m','Magnetic b.','Location','SouthEast')
%         legend('RE baricentro','Magnetic b.','Location','SouthEast')
%         %axis([t_i,t_f,min([PSI_magnetocentrum(:,2)' geo_bar]),max([PSI_magnetocentrum(:,2)' geo_bar])])
%         axis([t_start,t_end,0.8,max([PSI_magnetocentrum(:,2)' geo_bar])])
%         hold off
%         
%         
%         ax(4) = subplot(5,1,4);
%         hold on
%         temp = NEU213;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(temp.x(j1:end),temp.y(j1:end)./max(temp.y(j1:end)),'c','LineWidth',1.3);
%         grid on
%         hold on
%         temp = FBINA_M_RIS;
%         j_t0 = find(temp.x >= t0_offset,1);
%         offset = mean(temp.y(1:j_t0));
%         temp.y = temp.y - offset;
%         j1 = find(temp.x >= t_start,1);
%         j2 = find(temp.x >= t_end,1);
%         plot(mytime,temp.y(j1:j2),'k','LineWidth',1.5);
%         
%         set(gca,'FontSize',MyFontSize)
%         hold off
%         legend('neu213','Fbina','Location','SouthEast')
%         hold off
%         
%         
%         ax(5) = subplot(5,1,5);
%         plot(newtimepsi,PSI_magnetocentrum(:,3),'g');
%         grid on
%         set(gca,'FontSize',MyFontSize)
%         legend('PSI value','Location','SouthEast')
%         
%         linkaxes(ax,'x');
%         
%         figure(Nfigure+1)
%         plot(newtimepsi,(interp1(data.zs2.x,data.zs2.y,newtimepsi)+interp1(data.zs1.x,data.zs1.y,newtimepsi))/2.0,'b:','LineWidth',1.5);
%         grid on
%         set(gca,'FontSize',MyFontSize)
%         legend('Z baricentro +0.9m','Location','SouthEast')
%         
%         %pause
%         %close all
%         
%     end
%     
%     
% end

pause(mypause);

end
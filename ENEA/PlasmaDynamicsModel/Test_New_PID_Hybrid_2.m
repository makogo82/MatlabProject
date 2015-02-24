clc 
%close all
clear all

shot = 38869 %16 %38592 %38628 38628 27

Ts = 0.5E-3;
Tsampling = Ts;
t_start = 0.0;
t_end = 1.6;
data = CleanFTUdata_v6(t_start,t_end,Ts,shot,1);
mytime = data.time;

%%
%    theta = 0.65; 
%    plot(mytime,lsim(tf([(1-theta) 0],[1 -theta],Ts),data.dez/ ((1.5e-4) * (-1.0)),mytime),'b',...
%         mytime,data.dez/((1.5e-4) * (-1.0)),'r') 
%     grid on
%     
%     figure(6)
%     d_dez = mypseudo_derivative(mytime,data.dez,c_dez,d_dez);
%     plot(mytime,1000*data.dez,'r',mytime,data.Hcalc,'b',mytime,data.Hmis,'k',...
%          mytime,10*d_dez,'r--')
%     legend('dez','Hu','Hm','d_dez')
%     
    
%     %PID ENEA H
%     tau2 = 2.0E-3;
%     Kp_enea = 0.15%0.2199999988079;
%     Ki_enea = 1.0; %1.5
%     Kd_enea = 0.005;
%     
%     %xd(k) = (tau2*xd(k-1) + (u(k)-u(k-1)) )/(tau2+Ts)
%     %xd(z) = (z-1)/(z(Ts+tau2) -tau2) * u(z)
%     %xi(k) = xi(k-1) + 0.5*Ts*(u(k)+u(k-1))
%     %xi(z) = 0.5*Ts*(z+1)/(z-1)  *  u(z) 
%     PIDH_enea_z = (1/ (1.5e-4) * (-1.0))*(Kp_enea*tf([1],[1],Ts) + Kd_enea*tf([1 -1],[Ts+tau2, -tau2],Ts) +...
%                   Ki_enea*0.5*Ts*tf([1 1],[1 -1],Ts));
%     PIDH_enea_s = d2c(PIDH_enea_z,'zoh')
%     rlocus(PIDH_enea_s,1);
%     PIDH_enea_st = d2c(PIDH_enea_z,'tustin') 
%     hold on
%     rlocus(PIDH_enea_st,1);
%     legend('zoh','Tustin')
%     hold off
%     
%     PIDHsym = lsim(PIDH_enea_z,data.dez,mytime);
%     plot(mytime,data.Hcalc+17,'b',mytime,PIDHsym,'r')
%     legend('Hcalc','PID sym')
%     grid on
%     


%% Simulation
tau_filter = 0.75;
c_dez = 6;
d_dez = 18;
sat_ampl_pid = 350; %Ampere
sat_der_pid = 40.0; %Ampere/Ts
Kp =   0.18 %1.0* 5.0 *  0.22/Ts;
Ki = 1.0;
Kd =  0.0025;
Deadzone_der = 250;

Kd_dez  = 0*4E-5; %tra 6 e 4E-5
qKd_dez = 0*1;
Kd_2 = 0*5E-3;


switchingH = 0;
sigmaHpid1 = 1300.0; %600
sigmaHpid2 = -1000.0; %100
rampDerH = 4.0;
maxSatH = 200.0;     %130
errorDeazoneH = 20.0; %250
rampGainH = 0.00015;
rampGainHdez = 0.001;


c = c_dez;
d = d_dez;
y = data.dez/ ((1.5e-4) * (-1.0));
rif = 0.0*data.dez;
error = 0.0*data.dez;
d_error = zeros(length(data.dez),1);
mybuffer = zeros(d,1);
xfilter = zeros(length(data.dez),1);
xp = zeros(length(data.dez),1);
xd = zeros(length(data.dez),1);
xi = zeros(length(data.dez),1);
xpid = zeros(length(data.dez),1);
xpid_test = zeros(length(data.dez),1);
u_old = new_PID_HybridSat_1(1,error(1),tau_filter,Kp,Ki,Kd,c,d,sat_ampl_pid,sat_der_pid,Deadzone_der,Ts,sigmaHpid1,sigmaHpid2,rampDerH,maxSatH,errorDeazoneH,rampGainH,switchingH);


simula_mymodel = 1;

% A =[   1.01  -0.03138   0.02582
%        -0.02002    0.9728    0.3757
%        -0.02608  -0.07982   -0.2591];
% B = [-1.558e-05
%      -2.589e-05
%      -0.0001186];
% C = [ 0.8773  0.007708  -0.01387];
% D = 0;
%[A,B,C,D] = ssdata(c2d(ss(tf(-0.025,[9.07e-05,0.01215, 1.381 , 50])),Ts));
%[A,B,C,D] = ssdata(c2d(ss(A,B,C,D),Ts));

%mymodelDEZ = tf([0.001],[-1/15 1])*tf([1],[0.005 1]);
%[A,B,C,D] = ssdata(ss(c2d(mymodelDEZ,0.0005)));


%%

xplant = [10;-10];
My_disturbances = 0*(heaviside(mytime-0.5)-heaviside(mytime-0.55));

for j=2:length(mytime)
    
    %% MODELLO
    if(simula_mymodel==1)
        xplant = A*xplant + B*(u_old/1.5e-4 + My_disturbances(j));
        y(j) = C*xplant + D*u_old;
    end
    error(j) = y(j);
   
    %% CONTROLLORE
    
    xfilter(j) = tau_filter*xfilter(j-1) + (1-tau_filter)*error(j);
    
    if(j<d+1)
        mybuffer(j) = xfilter(j);
        d_error(j) = 0;
    else
        
        %PSEUDODERIVATIVE FUNCTION
        %----------------------------------
        %Buffer reset
        for k=1:d-1
            mybuffer(k) = mybuffer(k+1);
        end
        mybuffer(d) = xfilter(j);
        
        %ora calcolo le medie
        temp1 = 0;
        for k=1:c
            temp1 = temp1 + mybuffer(k);
        end
        temp1 = temp1/c;
        
        temp2 = 0;
        for k=d-c+1:d
            temp2 = temp2 + mybuffer(k);
        end
        temp2 = temp2/c;
        
        d_error(j) = (temp2-temp1) / (Ts*(d-c));
    end
    
    %proportional action
    xp(j) = Kp*xfilter(j);
    
    %derivative action with deadzone
    if(abs(d_error(j)) < Deadzone_der)
        xd(j) = 0;
    else
        xd(j) = Kd*d_error(j);
    end
    
    
    %partial pid, without the new value of the integrator
    xpid(j) = xp(j) + xd(j) + xi(j-1) ;
    
    %evaluation of amplitude and derivate saturations to obtain the integral new value
    temp1 = min( max(sat_ampl_pid - abs(xpid(j)),0.0),...
        max(sat_der_pid - abs(xpid(j)-xpid(j-1)),0.0));
    xi(j) = xi(j-1) + min(temp1,max(-temp1, Ki*(xfilter(j)+xfilter(j-1))*Ts/2.0));
    
    %final PID
    xpid(j) = min(sat_ampl_pid,max(-sat_ampl_pid, xpid(j-1)  + ...
        min(sat_der_pid,max(-sat_der_pid, -xpid(j-1) + xp(j) + xd(j) + xi(j))))) ;
    xpid_test(j) = new_PID_HybridSat_1(0,error(j),tau_filter,Kp,Ki,Kd,c,d,sat_ampl_pid,sat_der_pid,Deadzone_der,Ts,sigmaHpid1,sigmaHpid2,rampDerH,maxSatH,errorDeazoneH,rampGainH,switchingH,rampGainHdez,Kd_dez,qKd_dez,Kd_2);

    u_old = xpid_test(j);

end

if(simula_mymodel==1)
    
    figure(1)
    subplot(2,1,1)
    plot(mytime,y,'k');%,mytime,xpid_test,'r');
    ylim([-0.3,0.3])
    grid on
    legend('y')
    subplot(2,1,2)
    plot(mytime,xpid,'r',mytime,xpid_test,'b')
    grid on
    legend('xpid')
else
    
%     figure(1)
%     plot(mytime,y,'k',mytime,xfilter,'g',mytime,data.Hcalc+18,'r',mytime,xpid,'b',...
%         mytime,xpid_test,'b--')%,mytime,data.Hmis, 'm')
%     legend('y (dez)','xfilter','Hcalc','xpid','xpid_test')%,'Nmis')
%     grid on
%     axis([0,1.6,-150,150])
%    
%     
%     figure(2)
%     plot(mytime,3*data.dez,'k',mytime, xp,'m',mytime,xd,'r',mytime,xi,'g')
%     legend('dez','xp','xd','xi')
%     grid on
%     
    figure(3)
    ax(1) = subplot(3,1,1);
    plot(mytime,data.Hcalc,'r',mytime,xpid_test,'b',mytime,xpid,'g')
    legend('ALHExp','ALHRep','PIDRep');
    axis([0.1 1.5 -250 250]);
    title(num2str(shot));
   

    grid on
    ax(2) = subplot(3,1,2);
    plot(mytime,data.dez,'b',mytime,xfilter*((1.5e-4) * (-1.0)),'r', ...
        mytime,errorDeazoneH*((1.5e-4) * (-1.0)*ones(1,length(d_error))),'g--', ...
        mytime,-errorDeazoneH*((1.5e-4) * (-1.0)*ones(1,length(d_error))),'g--')
    axis([0.1 1.5  -0.2 0.2]);
    legend('DEZ','Filtered DEZ')
    
    grid on
    ax(3) = subplot(3,1,3);
    plot(mytime,-d_error,'b',mytime,sigmaHpid1*ones(1,length(d_error)),'k--',...
    mytime,sigmaHpid2*ones(1,length(d_error)),'g--',...
    mytime,-sigmaHpid2*ones(1,length(d_error)),'g--',...
    mytime,-sigmaHpid1*ones(1,length(d_error)),'k--');
    legend('d error');
    grid on
    %axis([0.1 1.5  1.03 1.28]);

    linkaxes(ax,'x');
    hold off
    
end
clc 
%close all
clear all
 
%% Shots for IDENTIFICATION
% FTU elongated    
%shots =  [39091    39092   39094 ...
%    38629 38628 38627 38616 38615 38614 38613 38612 38611 38610 38609 ...
%    38592 38591 38411 38410 38409 38404 38298 37930 37931 37932 ...
%    37870 37869];
 
shots =         [38591  39091   39092  39094  38629  38628 38627  38404  38616 38615];
shots_time = [    0.2   0.5     0.2    0.2    0.5    0.2   0.2    0.2    0.2   0.2 ];
shots_endtime = [ 1.29  0.904   1.15   1      0.99 1.195 1.22   1.12   0.857   1.6];
 
shot =38627
Ts = 0.5E-3;
samplingTime = Ts;
t_start = shots_time(find(shots==shot));
t_end = shots_endtime(find(shots==shot));
data = CleanFTUdata_v6(t_start,t_end,Ts,shot,1);
mytime = data.time;
%figure('Name','IPL')
%plot(mytime,data.IPLmis);
%% Simulation
tau_filter = 0.7; %0.75
c_dez = 2;
d_dez = 9;
c_ddez = 2;
d_ddez = 10;
 
satAmplPID = 350; %Ampere
satDerPID = 40.0; %Ampere/Ts
Kp = 0.18 %1.0* 5.0 *  0.22/Ts;
Ki = 1.0;
Kd = 0.0025;
derDeadzonePID = 50.0;
 
kd_error = 4E-5 % 4E-5 derivative term = (Kd + Kd_error*error )*d_error   
q_kd_error = 1.0 % to activate the Kd_error if e*d_e>0
kd_2 = 0.002    %0.002 piecewise gain added to Kd of the stardard derivative PID term when the |d_erro| is larger than SigmaHpid1
startTimeNonLinPID = 0.1;
endTimeNonLinPID = 1.5;
 
switchingH = 1; %1 se ramp-kick is ebabled
sigmaHpid1 = 5000; %activation threhold on error speed
sigmaHpid2 = -3000.0; %deactivation threhold on error speed
gammaHpid1 = 4.2E6; %activation threhold on error acceleration
gammaHpid2 = 10.0E5; %deactivation threhold on error acceleration
 
rampDerH = 4.0;
maxSatH = 300.0;     %130
errorDeazoneH = 90.0; %250
rampGainH = 0.00015
rampGainHdez = 0.001;
 
%%NEW PARAMETERS
acceleration_pid_H = 1; %if 1, than the hysteresis for the control is based on the acceleration instead of speed
avoid_chattering_PIDH = 10 % number of sampling time instants before the state can switch
avoid_chattering_PIDH_OFF = 5 % number of sampling time instants before the state can switch
theta_gain_0 = 0.2;
theta_gain_min =  0.0; %minimum value of the gain multiplier when switching-off the controller
theta_gain_max = 1.3; %maximum value of the gain multiplier when switching-off the controller
Kpe = 1.1* 0.20/150; % the nonlinear "proportional" gain Kp*(1+Kpe*error)
 
 
%% simulation vector data
y = data.dez/ ((1.5e-4) * (-1.0));
rif = 0.0*data.dez;
error_memo = 0.0*data.dez;
d_error_memo = zeros(length(data.dez),1);
dd_error_memo = zeros(length(data.dez),1);
xfilterState_memo = zeros(length(data.dez),1);
integralState_memo = zeros(length(data.dez),1);
xpidState_memo = zeros(length(data.dez),1);
xp_memo = zeros(length(data.dez),1);
xd_memo = zeros(length(data.dez),1);
qPIDH_memo = zeros(length(data.dez),1);
PIDsimple_memo  = zeros(length(data.dez),1);
switches_counter_memo  = zeros(length(data.dez),1);
switches_counter_OFF_memo  = zeros(length(data.dez),1);
qOFF_memo = zeros(length(data.dez),1);
theta_gains_memo = ones(length(data.dez),1); 
oscillationDezError_memo = zeros(length(data.dez),1);
switchCountDError_memo = zeros(length(data.dez),1);
switchCountDError=0;
 
 
%% MODELLO
modello_simulato_memo =  zeros(length(data.dez),2);
modello_amplificatore_memo = zeros(length(data.dez),2);
uu_memo =  zeros(length(data.dez),1);
number_of_resets_simulato = 1;%numero di volte si desidera far il reset del modello simulato
 
simula_mymodel = 1;
 
model_delay = 12; %sia INFERIORE A parti_modello, ritardo in numero di campioni
amplifier_positive_derivative_threshold = 25;
amplifier_negative_derivative_threshold = 15;
omega_p = 2*pi*100;
zita_p = 0.8;
a_2 = 0.2;
a_3 = 0.5E-5;
a_4  = 1.0E-6;
dis = 0.2;
s=tf([1 0],[1]);
[Ad,Bd,Cd,Dd] = ssdata(c2d(tf([omega_p^2],[1 2*zita_p*omega_p omega_p^2])*exp(-model_delay*s),Ts));
parti_modello = 2*model_delay;
 
%%
 
xplant = [10;-10];
My_disturbances = 0*(heaviside(mytime-0.5)-heaviside(mytime-0.55));
 
for j=2:length(mytime)
    
    %% MODELLO
    error_memo(j) = y(j);
    previousError = error_memo(j-1);
    
    %% Interfecciamento dati ingresso
    usecTime = data.time(j);
    xfilterState = xfilterState_memo(j-1);
    error = error_memo(j);
    qPIDH = qPIDH_memo(j-1);
    integralState = integralState_memo(j-1);
    xpidState = xpidState_memo(j-1);
    switches_counter = switches_counter_memo(j-1);
    switches_counter_OFF = switches_counter_OFF_memo(j-1);
    qOFF = qOFF_memo(j-1);
    theta_gains = theta_gains_memo(j-1);
    switchCountDError = switchCountDError_memo(j-1);
    oscillationDezError= oscillationDezError_memo(j-1);
    
    
    %% CODICE CONTROLLORE in Cpp
    
    if (usecTime >= startTimeNonLinPID &&  usecTime <= endTimeNonLinPID)
        NNLPID = 1;
    else
        NNLPID = 0;
    end
    
    xfilterState = tau_filter*xfilterState + (1-tau_filter)*error;
    d_error = on_line_mypseudoderivative(xfilterState,c_dez,d_dez,samplingTime);
    dd_error = on_line_mypseudoderivativeCopy(d_error,c_ddez,d_ddez,samplingTime);
    d_h = on_line_mypseudoderivativeH(data.Hmis(j),c_dez,d_dez,samplingTime);
    
    %% simulated model
    if(j > parti_modello)
        myinput = data.Hcalc(j-model_delay);
        %myinput = xpidState(j-model_delay);
        
        %modello_amplificatore_memo(j,1) = modello_amplificatore_memo(j-1,1) +  min(amplifier_derivative_threshold,max(-amplifier_derivative_threshold,...
        %    (Ad(1,:)*modello_amplificatore_memo(j-1,:)'+Bd(1)*myinput)-modello_amplificatore_memo(j-1,1)));
        %modello_amplificatore_memo(j,2) = modello_amplificatore_memo(j-1,2) +min(amplifier_derivative_threshold,max(-amplifier_derivative_threshold,...
        %    (Ad(2,:)*modello_amplificatore_memo(j-1,:)'+Bd(2)*myinput)-modello_amplificatore_memo(j-1,2) ));
        
        modello_amplificatore_memo(j,1) = modello_amplificatore_memo(j-1,1) +  (Ts*modello_amplificatore_memo(j-1,2));
        if(myinput-myinput_old<0)
            d_threshold = amplifier_negative_derivative_threshold;
        else
            d_threshold = amplifier_positive_derivative_threshold;
        end
        modello_amplificatore_memo(j,2) = min(d_threshold/Ts,max(-d_threshold/Ts, ...
            modello_amplificatore_memo(j-1,2) + Ts*(-2*zita_p*omega_p*modello_amplificatore_memo(j-1,2) - omega_p^2*modello_amplificatore_memo(j-1,1)+ omega_p^2*myinput) ));
        %uu_memo(j) = modello_amplificatore_memo(j,1);
        uu_memo(j) = data.Hcalc(j);
        myinput_old = myinput;
        modello_simulato_memo(j,1) = modello_simulato_memo(j-1,1)+ Ts*modello_simulato_memo(j-1,2);
        modello_simulato_memo(j,2) = modello_simulato_memo(j-1,2)+ Ts*(-a_2*modello_simulato_memo(j-1,2) + ...
            a_3*data.IPLmis(j)/abs(0.8-abs(modello_simulato_memo(j-1,1)))*uu_memo(j) + ...
            +sign(modello_simulato_memo(j-1,1))*(a_4*data.IPLmis(j)/(dis-min(dis*0.98,abs(modello_simulato_memo(j-1,1)))))^2);
    else
        myinput_old = 0;
        modello_amplificatore_memo(j,1) = data.Hmis(j);
        modello_amplificatore_memo(j,2) = d_h;%(data.Hmis(j)-data.Hmis(j-1))/Ts;
        modello_simulato_memo(j,1) = data.dez(j);
        modello_simulato_memo(j,2) = d_error*(((1.5e-4) * (-1.0)));%(data.dez(j)-data.dez(j-1))/Ts;
    end
    %%
    
    
    
    %proportional action
    xp = Kp*(1+Kpe*abs(xfilterState))*xfilterState;
    
    xd = 0.0;
    fabs_d_error = abs(d_error);
    fabs_xfilterState = abs(xfilterState);
       
    
    %derivative action with deadzone
    if(fabs_d_error > derDeadzonePID)
        %piecewise gain kd_2 that is added to Kd IF d_error is large
        kd_2temp = kd_2;
        if(fabs_d_error< sigmaHpid1) 
            kd_2temp = 0.0;
        end     
        if(d_error*xfilterState > 0.0)
            if(kd_2 > 0.0) 
                xd = (Kd + NNLPID*kd_2temp )* d_error; 
            else
                xd = (Kd + NNLPID*kd_error*fabs_xfilterState )* d_error;
            end
        else
            if(kd_2 > 0.0)
                xd = (Kd + NNLPID*(1.0-q_kd_error)*kd_2temp )* d_error;
            else
                xd = (Kd + NNLPID*(1.0-q_kd_error)*kd_error*fabs_xfilterState )* d_error;
            end
        end
    end
    
    
    
    %% Switching logic
    temp_qPIDH = qPIDH;
    temp_qOFF = qOFF;
    if(switchingH ~= 0)
       
        switch(qPIDH)
            
            case 0
                %Onset of kick is triggered by condition on error and
                %d_error
                
                
                if( ((xfilterState >= errorDeazoneH && abs(d_error)<= sigmaHpid2) ||...
                   (d_error*dd_error > 0 && (xp + xd + integralState)*d_error>0)) ...%the control is going to the opposite direction of the speed
                   && switches_counter_OFF >= avoid_chattering_PIDH_OFF ) 
                   qOFF = 0;
                    
                end
                 
                if(d_error >= sigmaHpid1 && xfilterState >= errorDeazoneH ...
                   && switches_counter >= avoid_chattering_PIDH )
                    qPIDH = 1;
                end
                if(d_error <= -sigmaHpid1 && xfilterState <= -errorDeazoneH ...
                   && switches_counter >= avoid_chattering_PIDH )
                    qPIDH = -1;
                end
                
                    
                
            case 1
                
                %setting the qPIDH variable
                if(d_error <= -sigmaHpid2  && switches_counter >= avoid_chattering_PIDH)
                    qPIDH = 0;
                end
                %setting qOFF to decide when "switch off" the Controller
                if(dd_error <= -gammaHpid2  && switches_counter_OFF >= avoid_chattering_PIDH_OFF)
                    qOFF = 1; %acceleration is "opposite" to speed, "switch-off" the controller
                end
                
                %if plasma continuous to go in the bad direction,
                %re-establish the old gains and skrink them less
                if(dd_error >= -gammaHpid2  && d_error >= sigmaHpid1 && switches_counter_OFF >= avoid_chattering_PIDH_OFF)
                    qOFF = 0;
                end
                
            case -1
                
                %setting the qPIDH variable
                if(d_error >= sigmaHpid2  && switches_counter >= avoid_chattering_PIDH)
                    qPIDH = 0;
                end
 
                                
                %setting qOFF to decide when "switch off" the Controller
                if(dd_error >= gammaHpid2  && switches_counter_OFF >= avoid_chattering_PIDH_OFF)
                    qOFF = 1; %acceleration is "opposite" to speed, "switch-off" the controller
                end
                %if(dd_error <= -gammaHpid1)
                %    
                %    qOFF = 0;
                %end
                %if plasma continuous to go in the bad direction,
                %re-establish the old gains and skrink them less
                if(dd_error <= gammaHpid2  && d_error <= -sigmaHpid1 && switches_counter_OFF >= avoid_chattering_PIDH_OFF)
                    qOFF = 0;
                end
                
        end
        
        
        %detect the inversion of qPIDH sign
        if ( qPIDH == 0 && temp_qPIDH ~= 0 && switchCountDError == 0)
            %begin counts and set qPIDHSwitch to qPIDH
            switchCountDError = switchCountDError + 1;
            qPIDHSwitch = temp_qPIDH;
        elseif( switchCountDError > 0)
            switchCountDError = switchCountDError +1;
            
            if (qPIDH == -qPIDHSwitch)
                switchCountDError = 0; %reset
                oscillationDezError = oscillationDezError + 1;
                %reduce the gain
            end
            if (qPIDH == qPIDHSwitch || switchCountDError > 100)
                switchCountDError = 0; %reset
                oscillationDezError = 0;
            end
        end
            
        
        
    else
        qPIDH = 0;
        qOFF = 0;
    end
    %avoid chattering
    if(temp_qPIDH ~= qPIDH)
         switches_counter = 0;    
    else
        switches_counter = min(avoid_chattering_PIDH, switches_counter + 1);
    end
    if(temp_qOFF ~= qOFF)
         switches_counter_OFF = 0;    
    else
        switches_counter_OFF = min(avoid_chattering_PIDH_OFF, switches_counter_OFF + 1);
    end
    
    %the gains of xp and xd is reduced by the factor theta_gains
    if(qOFF==1 )
        %the reducing factor is less if acceleration is large 
        theta_gains = 1.0*min(1.0, abs(dd_error)/(1.5*gammaHpid1)) + (1-min(1.0,abs(dd_error)/(1.5*gammaHpid1)))*theta_gain_0/(1+oscillationDezError);
    else
        theta_gains = 1.0;
    end
    
    
    %%
    if(qPIDH == 0)
        %partial pid, without the new value of the integrator
        xpid_new = xp + xd + integralState;
        %evaluation of amplitude and derivate saturations to obtain the integral new value
        temp1 = min( max(satAmplPID - abs(xpid_new),0.0),max(satDerPID - ...
               abs(xpid_new - xpidState),0.0));
        integralState = integralState + min(temp1,max(-temp1, ...
                          Ki*(xfilterState + previousError)*samplingTime*0.5));
        
        %final PID with amplitude and derivative saturation
        xpidState = min(satAmplPID,max(-satAmplPID, xpidState  + ...
            min(satDerPID,max(-satDerPID, -xpidState + xp + xd + integralState))));
        %xpidState(j) = min(satAmplPID,max(-satAmplPID, xpidState(j)  + ...
        %              min(satDerPID,max(-satDerPID, -xpidState(j) + xp(j) + xd(j) + integralState(j)))));    
    end
    
    if((qPIDH ~= 0 || (qOFF==1 && acceleration_pid_H)) && NNLPID > 0.0 )
        %adding the kick of the hybrid term
        
        if(acceleration_pid_H == 1)
            xpidState = min(maxSatH,max(-maxSatH, xpidState  + ...
            min(satDerPID,max(-satDerPID, -xpidState + theta_gains*(xp + xd) + integralState)))) ;
        else
            
            xpidState = min(maxSatH,max(-maxSatH, xpidState+ ...
                min(rampDerH,max(-rampDerH, qPIDH*rampGainH*fabs_d_error*(1+...
                rampGainHdez*fabs_xfilterState) )))) ;
        end
    end
    
    
    %% OUTPUT SIGNAL FOR SIMULATIONS
    xfilterState_memo(j) = xfilterState;
    error_memo(j) = error;
    d_error_memo(j) = d_error;
    dd_error_memo(j) = dd_error;
    xp_memo(j) = xp;
    xd_memo(j) = xd;
    qPIDH_memo(j) = qPIDH;
    integralState_memo(j) = integralState;
    xpidState_memo(j) = xpidState;
    switches_counter_memo(j) = switches_counter;
    switches_counter_OFF_memo(j) = switches_counter_OFF;
    qOFF_memo(j) = qOFF;
    theta_gains_memo(j) = theta_gains;
    oscillationDezError_memo(j)= oscillationDezError;
    switchCountDError_memo(j) = switchCountDError;
end
 
if(simula_mymodel==1)
    
    figure(1)
    ax(1) = subplot(2,1,1);
    plot(mytime,data.Hcalc,'b',mytime,data.Hmis,'r',mytime,modello_amplificatore_memo(:,1),'r--','LineWidth',2);%,mytime,xpid_test,'r');
    title(num2str(shot))
    grid on
    legend('ALHcalc','ALHmis','Estimated','Location','SouthWest');
    ax(2) = subplot(2,1,2);
    plot(mytime,data.dez,'r',mytime,modello_simulato_memo(:,1),'r--','LineWidth',2);
    grid on
    ylim([-0.2,0.2]);
    legend('data.dez','Estimated','Location','SouthWest');
    linkaxes(ax,'x');
    %%
    A = [0  1 0; 0 0 1; 0 0 0];
    rho = 1000;
    K = [3*rho; 3*rho^2; 1*rho^3];
    %eig(A-K*[1,0,0])
    XX = zeros(length(mytime),3);
    for i=2:length(mytime)
        XX(i,:) = XX(i-1,:) + Ts*(A*XX(i-1,:)' + K*(xfilterState_memo(i-1)-XX(i-1,1)))';
    end    
    %%
    figure(2)
    %c1 = a_3*data.IPLmis'./abs(0.8-abs(modello_simulato_memo(:,1))).*uu_memo;
    c_2 = 1E-5;
    c_3 = 0.8;
    b_pl = 0.29;% (data.zs1-data.zs2)/2.0;
    t1 = c_2*data.IPLmis'.*(1./(c_3-(data.dez+b_pl)) - 1./(-c_3-(data.dez+b_pl)) )'.*data.Hmis';
    c_4 = c_2*0.001;
    c_5 = 0.4;
    t2 = -c_4*data.Fmis'.*data.IPLmis'.*( 1./(c_5-(data.dez+b_pl)) + 1./(-c_5-(data.dez+b_pl)) )';
    %acc_int = lsim(tf([1],[1 0]),(t2+t1),mytime);
    
    ax(1) = subplot(2,1,1);
    plot(mytime, t1,'r',mytime, t2,'b',mytime,t1+t2,'g--',...%mytime,acc_int,'k--',...
        mytime,dd_error_memo*(((1.5e-4) * (-1.0))),'k',mytime,XX(:,3)*(((1.5e-4) * (-1.0))),'k:','LineWidth',2);
    grid on
    ylim([-3000, 3000]);
    title(num2str(shot))
    legend('c1','c2','c1+c2','dde','location','SouthWest');
    ax(2) = subplot(2,1,2);
    plot(mytime,d_error_memo,'lineWidth',2);
    grid on
    linkaxes(ax,'x');
    
    figure(3)
    ax(1) = subplot(4,1,1);
    plot(mytime,xfilterState_memo,'m-',mytime,XX(:,1),'m:',mytime,data.Hmis,'b',...
         mytime,data.Hcalc,'g--',mytime,xpidState_memo,'k--','LineWidth',2);
    title(num2str(shot))
    grid on
    legend('error','e estimate', 'Location','SouthWest');
    ax(2) = subplot(4,1,2);
    plot(mytime,d_error_memo,'b-',mytime,XX(:,2),'b:','LineWidth',2);
    grid on
    legend('d error','d estimate''Location','SouthWest');
    ax(3) = subplot(4,1,3);
    plot(mytime,dd_error_memo,'k-',mytime,XX(:,3),'k:','LineWidth',2);
    grid on
    legend('dd error','dd estimate','Location','SouthWest');
    ax(4) = subplot(4,1,4);
    plot(mytime,data.Fmis,'c-','LineWidth',2);
    grid on
    legend('I_F','Location','SouthWest');
    
    linkaxes(ax,'x');
   
else
    
    %%
    figure(1)
    ax(1) = subplot(4,1,1);
    plot(mytime,data.Hmis,'g',mytime,data.Hcalc,'k-',mytime,xfilterState_memo,'b-',mytime,xpidState_memo,'r-','LineWidth',2)
    legend('ALHmis','ALHcalc','erorr fit','PIDnew','Location','SouthWest');
    xlim([mytime(1),mytime(end)])% -250 250]);
    grid on
    title(num2str(shot));
   
    
    
    grid on
    ax(2) = subplot(4,1,2);
    plot(mytime,xfilterState_memo,'r',mytime,y,'r:', ...
        mytime,errorDeazoneH*ones(1,length(d_error_memo)),'g--', ...
        mytime,-errorDeazoneH*ones(1,length(d_error_memo)),'g--', 'LineWidth',2)
    %axis([0.1 1.5  -0.2 0.2]);
    legend('e filtered')   
    xlim([mytime(1),mytime(end)])
    grid on
    
    ax(3) = subplot(4,1,3);
    plot(mytime,d_error_memo,'b',mytime,qPIDH_memo*10000,'r',mytime,qOFF_memo*10000,'c-' , ...
        mytime,theta_gains_memo*10000,'m',mytime,sigmaHpid1*ones(1,length(d_error_memo)),'k--',...
    mytime,sigmaHpid2*ones(1,length(d_error_memo)),'g--',...
    mytime,-sigmaHpid2*ones(1,length(d_error_memo)),'g--',...
    mytime,-sigmaHpid1*ones(1,length(d_error_memo)),'k--','LineWidth',2);
    legend('d error','qPIDH','qOFF','qBouncing');
    grid on
    xlim([mytime(1),mytime(end)])
    %axis([0.1 1.5  1.03 1.28]);
     
    ax(4) = subplot(4,1,4);
    plot(mytime,dd_error_memo,'b',mytime,gammaHpid1*ones(1,length(d_error_memo)),'k--',...
    mytime,gammaHpid2*ones(1,length(d_error_memo)),'g--',...
    mytime,-gammaHpid2*ones(1,length(d_error_memo)),'g--',...
    mytime,-gammaHpid1*ones(1,length(d_error_memo)),'k--', 'LineWidth',2);
    legend('dd error');
    grid on
    xlim([mytime(1),mytime(end)])
    
    linkaxes(ax,'x');
    hold off
    %%
    figure(2)
    plot(mytime,xp_memo,'b:',mytime,xd_memo,'b--',mytime,integralState_memo,'b',...
         mytime,xpidState_memo,'r',mytime,xp_memo+xd_memo+integralState_memo,'k',...
         mytime,xfilterState_memo,'g',mytime,data.Hcalc,'c',mytime,data.Hmis,'c--','LineWidth',2)
    legend('xp','xd','xi','new','PID','error f','Hcalc','Hmis','Location','SouthWest');%,'FontSize',18)
    grid on
    title(num2str(shot));
    %%
    figure(3)
    theta0 = 1.2E-1;
    theta1 = 0.8;
    theta2 = 0;
    tau = 0.003;
    tau1 = 0.004;
    luky = theta2*error'-theta0*data.IPLmis./(theta1-abs(data.zs1+data.zs2)).*data.Hmis;
    plot(mytime,lsim(tf([1 0 0],[tau^2 2*tau 1]),error_memo,mytime),'b',mytime,luky,'r-',...
         mytime,lsim(tf([1],[ tau1 1]),luky,mytime),'k','LineWidth',2)
    legend('dd_error','luky','Location','SouthWest')
    grid on
    title(num2str(shot));
    %%
    figure(4)
    
    omega_z = 2*pi*60;
    zita = 0.2;
    s = tf([1 0],[1]);
    Ttau = 0.006;
    omega_p = 2*pi*120;
    zita_p = 0.7;
    Ttau_p = 0.00;
    plot(mytime,data.Hmis,'r-',mytime,data.Hcalc,'b-',...
         mytime,lsim(exp(-s*Ttau)*tf([1],[1/omega_z^2 2*zita/omega_z 1]),data.Hcalc,mytime),'k',...
         mytime,error_memo*0.4,'g',mytime,1.4*lsim(exp(-s*Ttau_p)*tf([1],[1/omega_p^2 2*zita_p/omega_p 1]),data.Hcalc,mytime),'g--', 'LineWidth',2)
    legend('Hmis','Hcalc','Hstimato','error','Location','SouthWest')
    grid on
    title(num2str(shot));
    %%
 
end
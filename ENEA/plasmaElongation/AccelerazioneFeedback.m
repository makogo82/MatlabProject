clc 
close all
clear all
 
%% Shots for IDENTIFICATION
% FTU elongated    
%shots =  [39091    39092   39094 ...
%    38629 38628 38627 38616 38615 38614 38613 38612 38611 38610 38609 ...
%    38592 38591 38411 38410 38409 38404 38298 37930 37931 37932 ...
%    37870 37869];
 
shots =         [ 38591  39091   39092  39094     38629       38628 38627  38404  38616 38615];
shots_time =    [ 0.2   0.7     1.0      0.85       0.85       0.89   1.10   0.2  0.752   0.2 ];
shots_endtime = [ 1.29  0.904   1.1        0.9      0.985       0.919  0.95   1.178   0.9  1.6];
 
shot =39092;
Ts = 0.5E-3;
samplingTime = Ts;
t_start = shots_time(find(shots==shot));
t_end = shots_endtime(find(shots==shot));
data = CleanFTUdata_v6(t_start,t_end,Ts,shot,1);
mytime = data.time;


%% Simulation
tau_filter = 0.7; %0.75
c_dez = 2;
d_dez = 9;
c_ddez = 2;
d_ddez = 10;
 
satAmplPID = 350; %Ampere
satDerPID = 40.0; %Ampere/Ts
Kp = 0.13; %1.0* 5.0 *  0.22/Ts;
Ki = 1.0;
Kd = 0.0025;
derDeadzonePID = 50.0;
 
kd_error = 4E-5 % 4E-5 derivative term = (Kd + Kd_error*error )*d_error   
q_kd_error = 1.0 % to activate the Kd_error if e*d_e>0
kd_2 = 0.002;    %0.002 piecewise gain added to Kd of the stardard derivative PID term when the |d_erro| is larger than SigmaHpid1
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
parti_modello = 2*model_delay;
 
%%
 
xplant = [10;-10];
My_disturbances = 0*(heaviside(mytime-0.5)-heaviside(mytime-0.55));
figure('Name','simulation')
plot(mytime,data.dez,'r','LineWidth',2); hold on;
d_error_dez = mypseudo_derivative(mytime, data.dez, c_ddez, d_ddez);
%%
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
    
    % AL-H model: y(t)=x(t-gamma); gamma=16 campioni (8ms)
    % Vertical Model: IPL and IF are constant
    % IH: is the input
    % Prediction of 100 samples (50 ms) 
    
    
    
    
    c_2H = 1.5E-5;
    c_2L = 0.5E-5;
    c_2 =  1E-5;
    c_3 = 0.7;
    c_4 = c_2*0.001;
    c_5 = 0.4
    Ts = 0.0005;
    
     
    IPLinput =  data.IPLmis(j);
    IFinput  = data.Fmis(j);
    
    myinput_old = 0;
    ampALH(1,1) = data.Hmis(j);
    ampALH(1,2) = d_h;%(data.Hmis(j)-data.Hmis(j-1))/Ts;
    

    if (j>2*model_delay) 
        xfilterStateP   = xfilterState_memo(j-1);
        integralStateP  = integralState_memo(j-1);
        xpidStateP = xpidState_memo(j-13:j-1)
        x1Ho(1) = data.dez(j);
        x2Ho(1) = d_error_dez(j);
        time(1) = usecTime;   
        HcalcInput(1:model_delay+1) = data.Hcalc(j-model_delay:j);
        xfilterStateP   = xfilterState_memo(j-1);
        xfilterStateStore   = xfilterState_memo(j-10:j-1);
        for jj = 1:10
            d_errorP = on_line_mypseudoderivativePrediction(xfilterStateStore(jj),c_dez,d_dez,samplingTime);
        end
        pid(1) = 0;
        for h = 1:1:100
            % begin prediction 
            
            %IHinput  = HcalcInput(h);  %get first element  
            %AMPLIFICATORE
            myinput =  HcalcInput(h);  %data.Hcalc(j-model_delay);
            ampALH(h+1,1) = ampALH(h,1) + (Ts*ampALH(h,2));
            if(myinput-myinput_old<0)
                d_threshold = amplifier_negative_derivative_threshold;
            else
                d_threshold = amplifier_positive_derivative_threshold;
            end
            
            ampALH(h+1,2) = min(d_threshold/Ts,max(-d_threshold/Ts, ...
            ampALH(h,2) + Ts*(-2*zita_p*omega_p*ampALH(h,2) - omega_p^2*ampALH(h,1)+ omega_p^2*myinput) ));
            myinput_old = myinput;
            IHinput = ampALH(h+1,1);
            
            
            [x1Ho(h+1), x2Ho(h+1)]  = verticalModel_Eulero(time(h), Ts, x1Ho(h), x2Ho(h), IPLinput, IHinput, IFinput, c_2H, c_3 , c_4 , c_5 );
            %                        
            errP =    x1Ho(h+1)  / ((1.5e-4) * (-1.0)); %% errore
            errPpre = x1Ho(h)    / ((1.5e-4) * (-1.0));
            
            xfilterStateP = tau_filter*xfilterStateP + (1-tau_filter)*errP;                        
            d_errorP = on_line_mypseudoderivativePrediction(xfilterStateP,c_dez,d_dez,samplingTime);
                       
            %proportional action
            xpP = Kp*(1+Kpe*abs(xfilterStateP))* xfilterStateP;
            xdP = Kd* d_errorP; 
            integralStateP =  integralStateP + Ki*(errP + errPpre)*samplingTime*0.5;
            pid(h+1) = xpP+xdP+integralStateP;          %PID
                        
            HcalcInput(model_delay+1+h) = pid(h+1);   % save input         

            time(h+1) = time(h)+Ts;

        end 
        
        bx(1) = subplot(2,1,1)        
        plot(time,x1Ho,'k:'); hold on;
        %ylim([-max(data.dez)*2 max(data.dez)*2]);
        grid on;
        bx(2) = subplot(2,1,2)
        %size([(time(1)-Ts*16):Ts:(time(1)-Ts) time])
   
        plot([(time(1)-Ts*model_delay):Ts:(time(1)-Ts) time],HcalcInput,'b--o'); hold on;
        
        data.Hcalc
        ylim([-500 500]);
        pause(0.01);
        grid on
        linkaxes(bx,'x')
        HcalcInput = [];
        x1Ho = []; 
        x2Ho = [];
        clear on_line_mypseudoderivativePrediction;
        %Vars=whos;
        %PersistentVars=Vars([Vars.global]);
        %PersistentVarNames={PersistentVars.name};
    end
      
    hold on;
    bx(1) = subplot(2,1,1)        
    plot(mytime,data.dez,'r'); hold on;
    %ylim([-max(data.dez)*2 max(data.dez)*2]);
    bx(2) = subplot(2,1,2)
    plot(mytime,data.Hcalc,'r'); hold on;
    ylim([-500 500])

    
    
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


figure('Name','If')
    plot(mytime,data.Fmis,'r'); hold on;

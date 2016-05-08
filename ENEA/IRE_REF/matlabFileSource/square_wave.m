function [deltaIpOut, enable]  = square_wave(time, activate, ipl_ref, ipl_ref_new, CQ_minIPL, plateauDetect,...
    pwm_amp, pwm_period, delay_activate_pwm, duty_cycle_pwm, enableBlock)

%square_wave(1.0, 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1);

%% COMPUTE PWM SIGNAL
% INPUT:
%           time     =
%           activate = boolean var for eneabled the alghoritm
%                      0 : active
%                      1 : non active
%           ipl      = absolute value of the  plasma current
%           dipl     = derivative of abs(ipl)
%           ipl_ref  = prep_reference
%
% OUTPUT:  y in R
%          y is the new reference for the ipl
%%

persistent    deltaIpl  startTime plateauDetectOld enableChase waitTime;

enable = 0;
deltaIpOut = 0;

if(activate && enableBlock)
    
    if(isempty(deltaIpl))
        deltaIpl = 0;
        waitTime = 0;
        enableChase = 0;
        plateauDetectOld= 0;
        startTime = 0;
    end
    
    
    % cq condition
    if(mod(plateauDetect,2) == 0 && plateauDetectOld ~= plateauDetect && plateauDetect >= 2)
        waitTime = time;
        enableChase = 1;
    end
    plateauDetectOld = plateauDetect;
    
    
    if((time-waitTime) >= delay_activate_pwm && enableChase >= 1 && plateauDetect >= 2)
        
        if(enableChase == 1)
            startTime =time;
        end        
        
        pwmTime = time-startTime;
        
        if( pwmTime < pwm_period )
            deltaIpl = -pwm_amp;
            enableChase = 2;
            if( pwmTime < pwm_period*duty_cycle_pwm )
                deltaIpl = pwm_amp;
                enableChase = 3;
            end
        else
            %reset pwmTime
            enableChase = 1;
        end   
    else
         deltaIpl = 0.0;          
    end
    

    if( pwm_amp > (ipl_ref_new-CQ_minIPL) || mod(plateauDetect,2) ~= 0)
        deltaIpl = 0;
        enableChase = 0;
    end
    
    
    enable = enableChase;
    deltaIpOut =  deltaIpl;
   
end


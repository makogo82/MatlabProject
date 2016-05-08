function [deltaIpOut, newipla_error, enable]  = ire_delta_ref(time, activate, vloop, ipl,ipl_sign, ipl_ref_new, plateauDetect,threshold_ipl_drift_apart, gamma_ipl_chase,delay_activate_ipl_chase, delay_activate_ipl_cq, threshold_vloop,threshold_ipl_vloop_activation,gamma_ipl_chase_vloop_mhd, enableBlock,CQ_minIPL)
%%ire_delta_ref(1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0,1.0, 1.0,1.0, 1.0, 1.0,1.0,1.0, 1.0,1.0)

%% COMPUTE NEW REFERENCE FOR RUNAWAY ELECTRON CURRENT BEAM WHEN ref and ipl are fall apart
% INPUT:
%           time     = 
%           activate = boolean var for eneabled the alghoritm
%                      0 : active
%                      1 : non active
%           ipl                         = absolute value of the  plasma current
%           dipl                        = derivative of abs(ipl)
%           ipl_ref                     = prep_reference
%           plateauDetect = 2,4,6...    =
%           threshold_ipl_fall_apart    =   
%           gamma_ipl_chase: 
%           delay_activate_ipl_chase:
%           enablePwm = 
% OUTPUT:  deltaIpOut: 
%          newipla_error: error 
%          enable
%%

persistent    deltaIpl   plateauDetectOld enableChase waitTime waitTimeChase;

enable = 0;
deltaIpOut = 0;
newipla_error = (ipl-ipl_ref_new);

if(activate && enableBlock)
    
    if(isempty(deltaIpl))
        deltaIpl = 0;
        waitTime = 0;
        enableChase = 0;
        plateauDetectOld= 0;
        waitTimeChase = 0.0;
    end
    
    if(plateauDetect >= 1)
        % cq condition
        if(mod(plateauDetect,2) == 0 && plateauDetectOld ~= plateauDetect)
            waitTime = time;
            enableChase = 1;
        end
        plateauDetectOld = plateauDetect;
        
        
        if(((time-waitTime) >= delay_activate_ipl_cq && 
             (time-waitTimeChase) >=delay_activate_ipl_chase) && 
              enableChase == 1 && ipl_ref_new >= CQ_minIPL)
            
            
            if((abs(newipla_error) >= threshold_ipl_drift_apart || -vloop*ipl_sign >= threshold_vloop) && mod(plateauDetect,2) == 0)
                enable = 1;
                newDelta = 0.0;
                
                if((abs(newipla_error) >= threshold_ipl_drift_apart))
                    newDelta = gamma_ipl_chase*newipla_error; 
                end
                
                %if vloop is due to mhd resistivity then 
                if( -vloop*ipl_sign >= threshold_vloop && newipla_error <= threshold_ipl_vloop_activation)
                    newDelta = -gamma_ipl_chase_vloop_mhd*ipl;
                end
                
                
                
                deltaIpl =  deltaIpl+newDelta;                
                waitTimeChase = time;
                
            end
            
            
            
        end
    end
    
    
    deltaIpOut =  deltaIpl;
    
end


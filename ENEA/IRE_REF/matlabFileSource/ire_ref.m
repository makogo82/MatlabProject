 function [Iref, enable, statusChangeNewIref] = ire_ref(time, activate, ipl, ipl_ref, plateauDetect, deltaTdes ,...
                                                        t_end ,beta, enableMultiplePlateau,...
                                                         percent_upgain_newiref,threshold_time_newiref,...
                                                         deltaIplChase, enableUpdateRefUp )
%% COMPUTE NEW REFERENCE FOR RUNAWAY ELECTRON CURRENT BEAM
% INPUT:
%           time     =
%           activate = boolean var for eneabled the alghoritm
%                      0 : active
%                      1 : non active
%           ipl      = absolute value of the  plasma current
%           dipl     = derivative of abs(ipl)
%           ipl_ref  = prep_reference
%           plateauDetect =
%           deltaTdes     =
%           t_end         = 
%           beta          =
%
% OUTPUT:  y in R
%          y is the new reference for the ipl
%%

persistent  enableCurrentDecay Iptpl tpl Ireftpl counter;

Iref = ipl_ref;
enable = 0;
statusChangeNewIref= 0;
if(activate)
    if(isempty(enableCurrentDecay))
        enableCurrentDecay = 0;
        Iptpl = 0;
        tpl = 0;
        Ireftpl = 0;
        counter = 0;
    end
    
    if(plateauDetect >=2)
        
        if( mod(plateauDetect,2) == 0 && ...
                enableCurrentDecay == (plateauDetect/2-1))
            % 2 && 0, 4 && 1 , 6 && 2 ecc..
            if(plateauDetect == 2 || enableMultiplePlateau == 1)
                %enableMultiplePlateau
                enableCurrentDecay = enableCurrentDecay + 1 ;
                %init var
                Iptpl   = ipl;
                if(enableCurrentDecay == 1)
                    %only ones
                    Ireftpl = ipl_ref;
                end
                tpl     = time;
                counter = 0;
            end
        end
        
        
        
        
    end
    
    
    
    if( enableCurrentDecay >= 1 )
        
        if(counter == 0)
            if(enableCurrentDecay == 1)
                Ireftpl = ipl_ref;
            end
        else
            
            %update program down
            newipla_error = (ipl-Ireftpl);
            
            %update program up
            if(enableUpdateRefUp == 1)
                if (newipla_error >= 0 && (t_end-time) >= threshold_time_newiref)
                    %update Iref_linear
                    Iptpl = ipl+(ipl*percent_upgain_newiref);
                    statusChangeNewIref=statusChangeNewIref + 1;
                end
            end
            
            
            
            % compute linear slop
            t_final = min(deltaTdes , t_end - tpl);
            alpha = Iptpl/(t_final);
            Iref_linear = (Iptpl) - alpha*(time-tpl) + deltaIplChase;
            % compute exp decay
            Ireftpl = (beta)*Ireftpl+(1-beta)*Iref_linear;
            
        end
        counter = counter +1;
        Iref =  Ireftpl;
        
        if(Iref <0) % disable
            enableCurrentDecay = -1;
        end
        
    end
    
    if(enableCurrentDecay == -1)
        Iref = 0;
    end
    

    
    enable = enableCurrentDecay;
    
end


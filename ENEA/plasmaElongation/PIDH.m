function [ xpidState , xp , xd , integralState] = PID_H( Kp, Kpe, Kd, error, previousError, d_error, xfilterState, integralState, xpidState, NNLPID, derDeadzonePID, satDerPID , satAmplPID, sigmaHpid1,kd_error,q_kd_error)

        %proportional action
        xfilterState = tau_filter*xfilterState + (1-tau_filter)*error;
        
        
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
  


end


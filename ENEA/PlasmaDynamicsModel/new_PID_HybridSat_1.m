function  PID_1 = new_PID_HybridSat_1(initall,error,tau_filter,Kp,Ki,Kd,c,d,sat_ampl_pid,sat_der_pid,Deadzone_der,Ts,sigmaHpid1,sigmaHpid2,rampDerH,maxSatH,errorDeazoneH,rampGainH,switchingH,rampGainHdez,kd_error,qKd_dez,kd_2)

%initall = 1 if the new_PID variables have to be initialized
%THE FIRST TIME THE FUNCTION IS CALL SET IT EGUAL TO 1
%c : the number of samples for avearaging (pseudoderivative)
%d : the number of samples distance  to evaluate the pseudoderivative
%sat_ampl_pid : amplitude saturation of the PID
%sat_der_pid : derivative amplitude saturation
%Kp : proportional gain
%Ki : integral gain
%Kd : derivative gain
%Deadzone_der : deadzone of the pseudoderivative
%Ts : sampling time


persistent old_error xi xpid mycounter mybuffer nullify xfilterState qPIDH kd_2_dead

if(initall==1)
    xfilterState = error;
    old_error = error;
    xi = 0.0;
    xpid = 0.0;
    mycounter = 0;
    mybuffer = zeros(d,1);
    qPIDH = 0;
    kd_2_dead = 0.0;
    if(c< 1 || d<2*c ||  d<2)
        disp('ERROR: c<1 or d < max(2*c,2) or  length(x) != length(y)')
        nullify = 0.0;
    else
        nullify = 1.0;
    end
    PID_1 = 0.0;
else
    
    
    mycounter = mycounter + 1;
    xfilterState = tau_filter*xfilterState + (1-tau_filter)*error;
    
    if(mycounter<d+1)
        
        mybuffer(mycounter) = xfilterState;
        d_error = 0;
        
    else
        %PSEUDODERIVATIVE FUNCTION
        %----------------------------------
        %Buffer reset
        for k=1:d-1
            mybuffer(k) = mybuffer(k+1);
        end
        mybuffer(d) = xfilterState;
        
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
        
        d_error = (temp2-temp1) / (Ts*(d-c));
    end
    
    
    
    
    if(switchingH == 1)
        %%setting the hybrid state jump (Switching logic)
        switch(qPIDH)
            case 0
                
                if(d_error >= sigmaHpid1 && xfilterState >= errorDeazoneH)
                    qPIDH = 1;
                end
                if(d_error <= -sigmaHpid1 && xfilterState <= -errorDeazoneH)
                    qPIDH = -1;
                end
                
            case 1
                if(d_error <= -sigmaHpid2  && d_error > -sigmaHpid1)
                    qPIDH = 0;
                end
                if(d_error <= -sigmaHpid1)
                    qPIDH = -1;
                end
            case -1
                if(d_error >= sigmaHpid2  && d_error < sigmaHpid1)
                    qPIDH = 0;
                end
                if(d_error >= sigmaHpid1)
                    qPIDH = 1;
                end
        end
        
    else
        qPIDH = 0;
    end
    
    %proportional action
    xp = Kp*xfilterState;
    
    %coefficiente della parte quadratica
    kd_2_dead = kd_2;
    if(abs(d_error) < sigmaHpid1)
        kd_2_dead = 0;
    end
    
    %derivative action with deadzone
    if(abs(d_error) < Deadzone_der)
        xd = 0;
    else
        if(xfilterState*d_error > 0.0)
          if(kd_2_dead > 0)
            xd = (Kd + kd_2_dead)*d_error; 
          else
            xd = (Kd + kd_error*abs(xfilterState))*d_error;
          end
        else
          if(kd_2_dead > 0)
            xd = (Kd + (1.0-qKd_dez)*kd_2_dead )*d_error;
          else
            xd = (Kd + (1.0-qKd_dez)*kd_error*abs(xfilterState) )*d_error;
          end
        end
        
    end
    
    
    if(qPIDH == 0 )
        
        %EVALUATE THE PID
        
        %partial pid, without the new value of the integrator
        xpid_new = xp + xd + xi ;
        
        %evaluation of amplitude and derivate saturations to obtain the integral new value
        temp1 = min( max(sat_ampl_pid - abs(xpid_new),0.0),...
            max(sat_der_pid - abs(xpid_new-xpid),0.0));
        xi = xi + min(temp1,max(-temp1, Ki*(xfilterState + old_error)*Ts/2.0));
        
        xpid = min(sat_ampl_pid,max(-sat_ampl_pid, xpid  + ...
            min(sat_der_pid,max(-sat_der_pid, -xpid + xp + xd + xi)))) ;
        %final PID
    end
    
    if(qPIDH ~= 0)
        xpid = min(maxSatH,max(-maxSatH, xpid  + ...
            min(rampDerH,max(-rampDerH, qPIDH*(rampGainH*abs(d_error)+rampGainHdez*abs(xfilterState)) )))) ;
        
    end
    
    
    
    old_error = xfilterState;
    derivativeState = d_error;
    PID_1 = xpid*nullify;
    
end

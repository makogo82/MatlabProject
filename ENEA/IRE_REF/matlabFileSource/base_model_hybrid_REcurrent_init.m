%% INIT 
Ts = 0.0002;                             % 200us
t0plasma  = 4.5;                         % zero plasma
t_end     = 2.0; 
rt2.system_period = 0.001;% discharge end
re.vloop_index = 119;
re.time_enable_ddint   = 0.1;
re.time_enable_new_ref = 0.1; 
re.time_activate = 0.1;
re.simtest       = 1;
%% IterativePseudoDerivative configuration
cq.deltaT = Ts;
cq.c = 2;
cq.d = 10;
cq.mediana = 0;
cq.reset = 0;
%% CQdetect
cq.minIPL = 40E3;   % always in abs value
cq.dIpl_threshold   = -1.1E7;
cq.counterThreshold = 3E-3/cq.deltaT; %3ms
%% new ref for ire current
beta      = 0.9;            % exponential decay rate
deltaTdes = 0.2;            % 300ms    
enableMultiplePlateau = 1;  % enable multiple plateau    
% shift the straight line of the new ipref
threshold_time_newiref=0.02;    %200ms
percent_upgain_newiref=0.1;    % if the (ipl current - newIreftpl) >= 0 then  increase Iref_linear of 5%
enableUpdateRefUp = 0;      %old version

%% chase of ipl_ref when abs(ipl-ipl_ref)>=threshold_ipl_drift_apart

trigger.threshold_ipl_time_fall_apart  = 0.1;
trigger.threshold_ipl_drift_apart      = 1e6;     %if this has to be really used the threshold_vloop MUST be very large 1e6 
trigger.gamma_ipl_chase                = 1;
trigger.enableTriggerChaseIpl          = 1;
trigger.delay_activate_ipl_chase       = 0.01;
trigger.delay_activate_ipl_cq          = 0.02;
trigger.threshold_vloop                = 2.5;     %E6; %if this has to be really used the threshold_ipl_drift_apart MUST be very large 1e6 
trigger.dipl_chase_threshold           = 1E7;
trigger.gamma_ipl_chase_vloop_mhd      = 0.1;   
trigger.threshold_ipl_vloop_activation = 2E4;     % if (ipl-ipl_ref) < threshold_ipl_vloop_activation then enable chase_vloop_mhd
%% pwm
pwm.enable_block   = 0;
pwm.amp            = 0.5E4;                % pwm amplitude
pwm.period         = 0.01;                 % 10ms
pwm.delay_activate = 0.005;                % delay of the pwm activation after plateaus detection
pwm.duty_cycle     = 0.5;                  % 50%


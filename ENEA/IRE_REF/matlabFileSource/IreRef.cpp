/*
    Author:mateusz.gospodarczyk@uniroma2.it
 */


#include "IreRef.h"
#include "LoadCDBObjectClass.h"

bool IreRef::ObjectLoadSetup(ConfigurationDataBase &cdbData, StreamInterface * err) {
    
    // Class Configuration
    
    if (!GCReferenceContainer::ObjectLoadSetup(cdbData, err)) {
        AssertErrorCondition(InitialisationError, "IreRef %s::ObjectLoadSetup: failed GCReferenceContainer constructor", Name());
        return False;
    }
    CDBExtended cdb(cdbData);
    samplingTime=0.0005;
    //--------------------------------------------- START CONF CQ DETECT BLOCK ------------------------------------------------------
    //CQIterativePseudoDerivative configuration
    /*
        cq.deltaT = Ts;
        cq.c = 2;
        cq.d = 10;
        cq.mediana = 0;
        cq.reset = 0;
     */
    pseudoDerivator = this->Find("CQ_IplDerivative");
    if (!pseudoDerivator.IsValid()) {
        AssertErrorCondition(InitialisationError, "IreRef %s::ObjectLoadSetup: do not find dIp error derivator ", Name());
        return False;
    }
    pseudoDerivator -> setSamplingTime(samplingTime);
    pseudoDerivator -> Reset();
    
    //CQdetect function
    
    
    cdb.ReadInt32(reverse_vloop_ipl, "Reverse_vloop_ipl",1); 																			//always in abs value
    AssertErrorCondition(Information, "IreRef %s::using reverse_vloop_ipl %f", Name(),reverse_vloop_ipl);
    
    
    cdb.ReadFloat(cq_minIPL,          "CQ_minIPL",40E3); 																			//always in abs value
    AssertErrorCondition(Information, "IreRef %s::using   minIPL %f", Name(), cq_minIPL);
    cdb.ReadFloat(cq_dIpl_threshold , "CQ_dIpl_threshold",-6E6);
    AssertErrorCondition(Information, "IreRef %s::using   cq_dIpl_threshold %f", Name(), cq_dIpl_threshold);
    cq_counterThreshold = (3.0E-3)/(samplingTime);
    
    cdb.ReadFloat(minIplActivationValue, "MinIplActivationValue",40E3); 																			//always in abs value
    AssertErrorCondition(Information, "IreRef %s::using   minIPL %f", Name(), minIplActivationValue);
    
    
    cdb.ReadFloat(activateTime, "ActivateTime",0.0); 																			//always in abs value
    AssertErrorCondition(Information, "IreRef %s::using activate %f", Name(), activateTime);
    
    //--------------------------------------------- END CONF CQ DETECT BLOCK ------------------------------------------------------
    
    
    
    //--------------------------------------------- START CONF PWM BLOCK ------------------------------------------------------
    cdb.ReadInt32(pwm_enable_block , 				"PWM_enable_block",0); //always in abs value
    cdb.ReadFloat(pwm_amp , 						"PWM_amp", 		0.5E4);												//pwm amplitude
    cdb.ReadFloat(pwm_period , 						"PWM_period"		, 0.01);										//pwm period 10 ms
    cdb.ReadFloat(pwm_delay_activate  , 			"PWM_delay_activate", 0.005);										//pwm delay of the pwm activation after plateaus detect
    cdb.ReadFloat(pwm_duty_cycle , 					"PWM_duty_cycle", 0.5);										//50%
    
    AssertErrorCondition(Information, "IreRef %s::using   pwm_enable_block %f", Name(), pwm_enable_block);
    AssertErrorCondition(Information, "IreRef %s::using   pwm_amp %f", Name(), pwm_amp);
    AssertErrorCondition(Information, "IreRef %s::using   pwm_period %f", Name(), pwm_period);
    AssertErrorCondition(Information, "IreRef %s::using   pwm_delay_activate  %f", Name(), pwm_delay_activate );
    AssertErrorCondition(Information, "IreRef %s::using   pwm_duty_cycle  %f", Name(), pwm_duty_cycle);
    //--------------------------------------------- END CONF PWM BLOCK ------------------------------------------------------
    
    
    
    
    //--------------------------------------------- START CONF  BLOCK ------------------------------------------------------
    
    // track of ipl_ref when abs(ipl-ipl_ref)>;//threshold_ipl_drift_apart
    cdb.ReadFloat( trigger_threshold_ipl_time_fall_apart  ,"TRIGGER_threshold_ipl_time_fall_apart"	,0.1);
    cdb.ReadFloat( trigger_threshold_ipl_drift_apart      ,"TRIGGER_threshold_ipl_drift_apart "   	,1e6);     //if this has to be really used the threshold_vloop MUST be very large 1e6
    cdb.ReadFloat( trigger_gamma_ipl_track                ,"TRIGGER_gamma_ipl_track"				,1);
    cdb.ReadInt32( trigger_enable_trigger_track_ipl       ,"TRIGGER_enable_trigger_track_ipl"       ,1);
    cdb.ReadFloat( trigger_delay_activate_ipl_track       ,"TRIGGER_delay_activate_ipl_track"		,0.01);
    cdb.ReadFloat( trigger_delay_activate_ipl_cq          ,"TRIGGER_delay_activate_ipl_cq"			,0.02);
    cdb.ReadFloat( trigger_threshold_vloop                ,"TRIGGER_threshold_vloop"				,2.5);     //E6; //if this has to be really used the threshold_ipl_drift_apart MUST be very large 1e6
    cdb.ReadFloat( trigger_dipl_track_threshold           ,"TRIGGER_dipl_track_threshold"			,1E7);
    cdb.ReadFloat( trigger_gamma_ipl_track_vloop_mhd      ,"TRIGGER_gamma_ipl_track_vloop_mhd" 		,1);
    cdb.ReadFloat( trigger_threshold_ipl_vloop_activation ,"TRIGGER_threshold_ipl_vloop_activation"	,2E4);     // if (ipl-ipl_ref) < threshold_ipl_vloop_activation then enable track_vloop_mhd
    
    
    
    //--------------------------------------------- START CONF  BLOCK ------------------------------------------------------
    //new ref for ire current
    cdb.ReadFloat( 	ireRef_beta      				,"IREREF_beta", 						0.9);//            % exponential decay rate
    cdb.ReadFloat( 	ireRef_deltaTdes 				,"IREREF_deltaTdes", 				0.2);//            % 300ms
    cdb.ReadInt32( 	ireRef_enableMultiplePlateau 	,"IREREF_enableMultiplePlateau", 	1);//  % enable multiple plateau
    //shift the straight line of the new ipref
    cdb.ReadFloat( 	ireRef_threshold_time_newiref 	,"IREREF_threshold_time_newiref", 	0.02);//    %200ms
    cdb.ReadFloat( 	ireRef_percent_upgain_newiref 	,"IREREF_percent_upgain_newiref",	 0.1);//    % if the (ipl current - newIreftpl) >= 0 then  increase Iref_linear of 5%
    cdb.ReadFloat( 	ireRef_enableUpdateRefUp 		,"IREREF_enableUpdateRefUp",			 0);//      %old version
    
    
    cdb.ReadFloat( 	t0plasma 	,"T_zero_plasma",	 10);
    cdb.ReadFloat( 	t_end 		,"T_end",            12);
    
    bool res = True;
    
    return res;
}


//main function that compute the new reference
float IreRef::new_ire_ref(int32 time, float ipl, float ipl_ref, float vloop, ControllerGAMOutputStructure* outputData){
    // execution order
    // 1:enable cq and palteau detection
    outputData->ipl_ref_old = ipl_ref;
    if(cq_detect_output != -1)
        ipl_sign = (ipl > 0.0) ? 1 : ((ipl < 0.0) ? -1 : 0); //sign
    
    timeRel = (float)(time*1E-6-t0plasma);    
    ipl = fabs(ipl);
    ipl_ref = fabs(ipl_ref);
    
    outputData->time = timeRel;
    
    if(time >= activateTime){

        activate = 1;
        
        dIp = pseudoDerivator -> Evaluate(timeRel, ipl);
        outputData->dIp = dIp;
            
            
            outputData->ipl_sign = ipl_sign;
            
            
            cq_detect(ipl);
            outputData->cq_detect_output = cq_detect_output;     //ok
            
            if(switch_new_iref){
                ipl_ref = iref_output;
            }
            
            
            // 2 compute delta_ipla by ire_delta_ref function
            ire_delta_ref( timeRel, vloop,  ipl, ipl_ref);          // delta ipl compute
            outputData->enable_delta_ref = enable_delta_ref;
            // 3
            ire_ref(timeRel, ipl, ipl_ref);            
            outputData->enable_current_tracking =enable_current_tracking;            

            
            //deltaIpOut
            // 4 conmpute pwm if the function is enabled
            square_wave( timeRel,  ipl_ref);
            outputData->enable_pwm = enable_square;            
            outputData->ipl_error = newipla_error;

            
            //diagnostic only
            outputData->iref_output    = iref_output*ipl_sign;
            outputData->deltaIpOut     = deltaIpOut;
            outputData->deltaIpOut_pwm = (deltaIpOut_pwm + iref_output) * ipl_sign;
            outputData->iref_linear    = deltaIpOut_pwm;
            
            //recover the right sign of the signal
            return iref_output*ipl_sign;
            
            
        
    }
    
    return  ipl_ref;
    
    
}




// new algorithm for current quench and plateau detection
void IreRef::cq_detect(float ipl){
    
    // % CurrentQuench and RE palteau detection
    //  INPUT:
    //           activate = boolean var for eneabled the alghoritm
    //                       0 : active
    //                       1 : non active
    //            ipl = absolute value of the  plasma current
    //            dipl, minIPL, dipl_threshold, counterThreshold : min
    //
    //  OUTPUT:  y in R
    //          -1 : detection non active
    //           0 : detection active
    //           1 : CQ detect
    //           2 : RE palteau detect
    // %
    // Arguments    : int activate
    //                float ipl
    //                float dipl
    //                float minIPL
    //                float dipl_threshold
    //                float counterThreshold
    // Return Type  : float
    
    
    if (activate == 1 ) {
        // STEP1: FIND CQ
        if (ipl >= cq_minIPL) {
            
            if ((dIp <= cq_dIpl_threshold) && fmod(startDetection, 2.0) == 0.0)
                startDetection++;
            
            
            if ((dIp >= cq_dIpl_threshold) && fmod(startDetection, 2.0)  == 1.0)
                startDetection++;
            
        } else {
            if (startDetection > 0.0) {
                //finalLoss
                startDetection = -1.0;
            }
        }
        
        cq_detect_output = startDetection;
    } else {
        cq_detect_output = -2;
    }
    
    
    //cq_detect_output is the output
    
}

//square wave for pwm
void IreRef::square_wave(float time, float  ipl_ref_new){
    
    // % COMPUTE PWM SIGNAL
    //  INPUT:
    //            time     =
    //            activate = boolean var for eneabled the alghoritm
    //                       0 : active
    //                       1 : non active
    //            ipl      = absolute value of the  plasma current
    //            dipl     = derivative of abs(ipl)
    //            ipl_ref  = prep_reference
    //
    //  OUTPUT:  y in R
    //           y is the new reference for the ipl
    // %
    

    deltaIpOut_pwm = 0.0;

    if ((activate == 1 ) && (pwm_enable_block == 1)) {
        
        //  cq condition
        if ( fmod(cq_detect_output, 2) == 0 && (plateauDetectOld_pwm != cq_detect_output) && (cq_detect_output >= 2) ) {
            waitTime_pwm = time;            
            enable_square = 1;
        }
        
        plateauDetectOld_pwm = cq_detect_output;
        
        if ((time - waitTime_pwm >= pwm_delay_activate) && (enable_square >= 1) && (cq_detect_output >= 2)) {
            
            if (enable_square == 1) {
                startTime_pwm = time;
            }
            
            pwmTime = time - startTime_pwm;
            
            if (pwmTime < pwm_period) {
                deltaIpOut_pwm= -pwm_amp;
                enable_square= 2;
                if (pwmTime < pwm_period * pwm_duty_cycle) {
                    deltaIpOut_pwm = pwm_amp;
                    enable_square= 3;
                }
            } else {
                // reset pwmTime
                enable_square= 1;
            }
        } else {
            deltaIpOut_pwm = 0.0;                   
        }
        
        if ((pwm_amp > ipl_ref_new - cq_minIPL) || fmod(cq_detect_output , 2)  != 0) {
            deltaIpOut_pwm = 0.0;
            enable_square = -1;         //stop
        }
        
    }
     

    
}

// ire delta_ref
void IreRef::ire_delta_ref(float time, float vloop, float ipl, float ipl_ref_new){
    
    // use 	cq_detect_output
    // % COMPUTE NEW REFERENCE FOR RUNAWAY ELECTRON CURRENT BEAM WHEN ref and ipl are fall apart
    //  INPUT:
    //            time     =
    //            activate = boolean var for eneabled the alghoritm
    //                       0 : active
    //                       1 : non active
    //            ipl                         = absolute value of the  plasma current
    //            dipl                        = derivative of abs(ipl)
    //            ipl_ref                     = prep_reference
    //            plateauDetect = 2,4,6...    =
    //            threshold_ipl_fall_apart    =
    //            gamma_ipl_track:
    //            delay_activate_ipl_track:
    //            enablePwm =
    //  OUTPUT:  deltaIpOut:
    //           newipla_error: error
    //           enable
    // %
    
    
    // OUTPUT VAR
    enable_delta_ref = 0;
    deltaIpOut = 0.0;
    newipla_error = ipl - ipl_ref_new;
    
    
    if ((activate ==1 ) && (trigger_enable_trigger_track_ipl == 1) && ipl >= minIplActivationValue && (cq_detect_output >= 2) ) {
        if (cq_detect_output  >= 1) {
            //  cq condition
            if (fmod(cq_detect_output,2) == 0 &&  (plateauDetectOld != cq_detect_output )) {
                waitTime = time;
                enable_tracking = 1;
                enable_delta_ref = 1.5;
            }
            
            plateauDetectOld = cq_detect_output ;
            
            
            if ((time - waitTime >= trigger_delay_activate_ipl_cq) && (time - waitTimetrack >=
            trigger_delay_activate_ipl_track) && (enable_tracking == 1) && (ipl_ref_new >= cq_minIPL)) {
                enable_delta_ref = 2;
                
                if (fabs(newipla_error) >= trigger_threshold_ipl_drift_apart || (reverse_vloop_ipl * vloop *ipl_sign) >= trigger_threshold_vloop && fmod(cq_detect_output,2) == 0) {
                    
                    enable_delta_ref = 3;
                    newDelta = 0.0;
                    if (fabs(newipla_error) >= trigger_threshold_ipl_drift_apart) {
                        newDelta = trigger_delay_activate_ipl_track * newipla_error;                        
                        enable_delta_ref = 4;
                    }
                    
                    // if vloop is due to mhd resistivity then
                    if(reverse_vloop_ipl* vloop * ipl_sign >= trigger_threshold_vloop && newipla_error <= trigger_threshold_ipl_vloop_activation){
                        //newDelta = -trigger_gamma_ipl_track_vloop_mhd * ipl;                            
                        newDelta = -trigger_gamma_ipl_track_vloop_mhd*fabs(newipla_error); 
                        enable_delta_ref = 5;
                    }
                    
                    
                    deltaIpl += newDelta;
                    waitTimetrack = time;
                    
                    
                }
            }
            
            deltaIpOut = deltaIpl;
            
        }
        
    }
}
//ire_ref
void IreRef::ire_ref(float time, float ipl, float ipl_ref){
    
  /*
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
   */
    
    
    
    iref_output = ipl_ref;
    statusChangeNewIref = 0.0;
    enable_current_tracking = 0;
    
    if (activate == 1) {
        if ((cq_detect_output >= 2) || (cq_detect_output == -1)) {
            
            
            
            
            switch_new_iref = true;
            if( fmod(cq_detect_output, 2) == 0 && current_tracking_status == (cq_detect_output/2-1)){
                
                if (((cq_detect_output == 2) || (ireRef_enableMultiplePlateau == 1))) {
                    
                    //  2 && 0, 4 && 1 , 6 && 2 ecc..
                    // enableMultiplePlateau
                    current_tracking_status++;
                    // init var
                    Iptpl = ipl;
                    if (current_tracking_status == 1) {
                        // only ones
                        Ireftpl = ipl_ref;
                    }
                    
                    tpl = time;
                    counter = 0.0;
                }
            }
            
            
        }
        
        
        
        
        if (current_tracking_status >= 1) {
            if (counter == 0) {
                if (current_tracking_status == 1) {
                    Ireftpl = ipl_ref;
                }
            } else {
                // update program down
                // update program up
                
                t_final = MIN(ireRef_deltaTdes, t_end - tpl);
                
                ireRef_alpha = Iptpl/(t_final);
                Iref_linear = (Iptpl) - ireRef_alpha*(time-tpl) + deltaIpOut; // deltaIpOut from previouse function
                //compute exp decay
                
                if(cq_detect_output == -1){ //final loss
                    Iref_linear = 0;                    
                }
                Ireftpl = (ireRef_beta)*Ireftpl+(1-ireRef_beta)*Iref_linear;
                
                
                   
                
                
            }
            
            counter++;
            
            iref_output = Ireftpl;
            
            if (Ireftpl < 0.0) {
                //  disable
                current_tracking_status = -1;
            }
        }
        
        if (current_tracking_status == -1) {
            iref_output = 0.0;
        }
        enable_current_tracking= current_tracking_status;
    }
    
    
}





OBJECTLOADREGISTER(IreRef, "$Id:IreRef.cpp,v 1.1.1.1 2010-01-20 12:26:47 gangapc Exp $")


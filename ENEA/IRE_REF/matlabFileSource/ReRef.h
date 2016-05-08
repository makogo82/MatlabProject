/*
    Author:mateusz.gospodarczyk@uniroma2.it
 */

#ifndef IreRefOutput_H_
#define IreRefOutput_H_

//#include "GAM.h"
#include "GCRTemplate.h"
#include "GCReferenceContainer.h"
#include "PseudoDerivator.h"
#include "ControllerGAMOutputStructure.h"
#include <math.h>

inline float MAX(float a, float b) {
    return a > b ? a : b;
}

inline float MIN(float a, float b) {
    return a < b ? a : b;
}



OBJECT_DLL(IreRef)
class IreRef : public GCReferenceContainer {
    OBJECT_DLL_STUFF(IreRef)
    // Parameters
    private:
        
        
        GCRTemplate<PseudoDerivator> pseudoDerivator;
        
        float activateTime;    
    	float t0plasma;
    	float t_end 	;
        float timeRel;
       
        int32 ipl_sign;
        int32 reverse_vloop_ipl;
        bool switch_new_iref;

        float t_final;
        float dIp;
        float minIplActivationValue;
        int32 activate;
        float samplingTime;
        
        float cq_minIPL;// 40E3;   // always in abs value
        float cq_dIpl_threshold   ;// -1.1E7;
        float cq_counterThreshold ;// 3E-3/cq.deltaT; //3ms
        int32 cq_detect_output;
        
        
        // pwm parameter
        int32 pwm_enable_block;//   ;// 0;
        float pwm_amp;//            ;// 0.5E4;                // pwm amplitude
        float pwm_period;//         ;// 0.01;                 // 10ms
        float pwm_delay_activate;// ;// 0.005;                // delay of the pwm activation after plateaus detection
        float pwm_duty_cycle;//     ;// 0.5;                  // 50//
        float pwmTime;
        //output
        int enable_pwm;
        float deltaIpOut_pwm;
        
        // track of ipl_ref when abs(ipl-ipl_ref)>;//threshold_ipl_drift_apart
        float trigger_threshold_ipl_time_fall_apart  ;// 0.1;
        float trigger_threshold_ipl_drift_apart      ;// 1e6;     //if this has to be really used the threshold_vloop MUST be very large 1e6
        float trigger_gamma_ipl_track                ;// 1;
        int32 trigger_enable_trigger_track_ipl    ;// 1;
        float trigger_delay_activate_ipl_track       ;// 0.01;
        float trigger_delay_activate_ipl_cq          ;// 0.02;
        float trigger_threshold_vloop                ;// 2.5;     //E6; //if this has to be really used the threshold_ipl_drift_apart MUST be very large 1e6
        float trigger_dipl_track_threshold           ;// 1E7;
        float trigger_gamma_ipl_track_vloop_mhd      ;// 0.1;
        float trigger_threshold_ipl_vloop_activation ;// 2E4;     // if (ipl-ipl_ref) < threshold_ipl_vloop_activation then enable track_vloop_mhd
        float newipla_error;
        float deltaFinalLoss;
        
        //new ref for ire current
        float ireRef_beta;						// = 0.9;     % exponential decay rate
        float ireRef_alpha;						// = 0.9;     % exponential decay rate
        float ireRef_deltaTdes;				// = 0.2;     % 300ms
        int32 ireRef_enableMultiplePlateau;	// = 1;  	  % enable multiple plateau
        //shift the straight line of the new ipref
        float ireRef_threshold_time_newiref;	//= 0.02;     %200ms
        float ireRef_percent_upgain_newiref;	//= 0.1;   	  % if the (ipl current - newIreftpl) >= 0 then  increase Iref_linear of 5%
        float ireRef_enableUpdateRefUp;		//= 0.0;      %old version
        float iref_output;
        
        //
        float deltaIpl;
        float waitTime_pwm;
        int32 enable_square;
        int32 plateauDetectOld_pwm;
        float startTime_pwm;
        float newDelta;
        int enable_delta_ref;
        int enable_current_tracking;
        //
        float startDetection;
        //
        int32 current_tracking_status;
        
        float Iptpl;
        float tpl;
        float Ireftpl;
        int32 counter;
        int enable_tracking;
        //
        
        float waitTime;
        float enabletrack;
        int32 plateauDetectOld;
        float startTime;
        float waitTimetrack;
        
        float statusChangeNewIref;
        float enable_ire_ref;
        float deltaIpOut;
        float Iref_linear;
        
        
        public:
            
            IreRef() {
                switch_new_iref = false;
                square_wave_init();
                ire_ref_init();
                ire_delta_ref_init();
                cq_detect_init();
                
            }
            
            void square_wave_init()
{
    
    waitTime_pwm = 0.0;
    enable_square = 0;
    plateauDetectOld_pwm = 0;
    startTime_pwm = 0.0;
            }
            
            void ire_ref_init()
{
    current_tracking_status = 0;
    Iptpl = 0.0;
    tpl = 0.0;
    Ireftpl = 0.0;
    counter = 0;
            }
            
            void ire_delta_ref_init()
{
    
    enable_tracking = 0;
    deltaIpl = 0.0;
    waitTime = 0.0;
    enabletrack = 0.0;
    plateauDetectOld = 0;
    waitTimetrack = 0.0;
            }
            
            void cq_detect_init()
{
    startDetection = 0.0;
            }
            
            ~IreRef() {
            }
            
            // new algorithm for current quench and plateau detection
            void cq_detect(float ipl);
            
            //square wave for pwm
            void square_wave(float timeRel, float  ipl_ref_new);
            
            // ire delta_ref
            void ire_delta_ref(float timeRel, float vloop, float ipl, float ipl_ref_new);
            
            //ire_ref
            void ire_ref(float timeRel, float ipl, float ipl_ref);
            
            //main function
            float new_re_ref(int32 time, float ipl, float ipl_ref, float vloop, ControllerGAMOutputStructure* outputData);
            
            virtual bool ObjectLoadSetup(ConfigurationDataBase &cdbData, StreamInterface * err);
            
};
#endif

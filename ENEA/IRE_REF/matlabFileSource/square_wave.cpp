//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: square_wave.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 16:17:07
//

// Include Files
#include "rt_nonfinite.h"
#include "square_wave.h"

// Variable Definitions
static double deltaIpl;
static double startTime;
static double plateauDetectOld;
static double enableChase;
static double waitTime;

// Function Definitions

//
// square_wave(1.0, 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1);
// Arguments    : double time
//                int activate
//                double ipl_ref
//                double ipl_ref_new
//                double CQ_minIPL
//                double plateauDetect
//                double pwm_amp
//                double pwm_period
//                double delay_activate_pwm
//                double duty_cycle_pwm
//                int enableBlock
//                double *deltaIpOut
//                double *enable
// Return Type  : void
//
void square_wave(double time, int activate, double, double ipl_ref_new, double
                 CQ_minIPL, double plateauDetect, double pwm_amp, double
                 pwm_period, double delay_activate_pwm, double duty_cycle_pwm,
                 int enableBlock, double *deltaIpOut, double *enable)
{
  double pwmTime;

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
  *enable = 0.0;
  *deltaIpOut = 0.0;
  if ((activate != 0) && (enableBlock != 0)) {
    //  cq condition
    if ((plateauDetect - floor(plateauDetect / 2.0) * 2.0 == 0.0) &&
        (plateauDetectOld != plateauDetect) && (plateauDetect >= 2.0)) {
      waitTime = time;
      enableChase = 1.0;
    }

    plateauDetectOld = plateauDetect;
    if ((time - waitTime >= delay_activate_pwm) && (enableChase >= 1.0) &&
        (plateauDetect >= 2.0)) {
      if (enableChase == 1.0) {
        startTime = time;
      }

      pwmTime = time - startTime;
      if (pwmTime < pwm_period) {
        deltaIpl = -pwm_amp;
        enableChase = 2.0;
        if (pwmTime < pwm_period * duty_cycle_pwm) {
          deltaIpl = pwm_amp;
          enableChase = 3.0;
        }
      } else {
        // reset pwmTime
        enableChase = 1.0;
      }
    } else {
      deltaIpl = 0.0;
    }

    if ((pwm_amp > ipl_ref_new - CQ_minIPL) || (plateauDetect - floor
         (plateauDetect / 2.0) * 2.0 != 0.0)) {
      deltaIpl = 0.0;
      enableChase = 0.0;
    }

    *enable = enableChase;
    *deltaIpOut = deltaIpl;
  }
}

//
// Arguments    : void
// Return Type  : void
//
void square_wave_init()
{
  deltaIpl = 0.0;
  waitTime = 0.0;
  enableChase = 0.0;
  plateauDetectOld = 0.0;
  startTime = 0.0;
}

//
// File trailer for square_wave.cpp
//
// [EOF]
//

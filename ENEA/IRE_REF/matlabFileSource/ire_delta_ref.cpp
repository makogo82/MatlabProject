//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ire_delta_ref.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 16:09:20
//

// Include Files
#include "rt_nonfinite.h"
#include "ire_delta_ref.h"

// Variable Definitions
static double deltaIpl;
static double plateauDetectOld;
static double enableChase;
static double waitTime;
static double waitTimeChase;

// Function Declarations
static double rt_roundd_snf(double u);

// Function Definitions

//
// Arguments    : double u
// Return Type  : double
//
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

//
// %ire_delta_ref(1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0,1.0, 1.0,1.0, 1.0, 1.0,1.0,1.0, 1.0,1.0)
// Arguments    : double time
//                int activate
//                double vloop
//                double ipl
//                int ipl_sign
//                double ipl_ref_new
//                double plateauDetect
//                double threshold_ipl_drift_apart
//                double gamma_ipl_chase
//                int delay_activate_ipl_chase
//                double delay_activate_ipl_cq
//                double threshold_vloop
//                double threshold_ipl_vloop_activation
//                double gamma_ipl_chase_vloop_mhd
//                int enableBlock
//                double CQ_minIPL
//                double *deltaIpOut
//                double *newipla_error
//                double *enable
// Return Type  : void
//
void ire_delta_ref(double time, int activate, double vloop, double ipl, int
                   ipl_sign, double ipl_ref_new, double plateauDetect, double
                   threshold_ipl_drift_apart, double gamma_ipl_chase, int
                   delay_activate_ipl_chase, double delay_activate_ipl_cq,
                   double threshold_vloop, double threshold_ipl_vloop_activation,
                   double gamma_ipl_chase_vloop_mhd, int enableBlock, double
                   CQ_minIPL, double *deltaIpOut, double *newipla_error, double *
                   enable)
{
  boolean_T guard1 = false;
  double d0;
  int i0;
  double newDelta;

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
  //            gamma_ipl_chase:
  //            delay_activate_ipl_chase:
  //            enablePwm =
  //  OUTPUT:  deltaIpOut:
  //           newipla_error: error
  //           enable
  // %
  *enable = 0.0;
  *deltaIpOut = 0.0;
  *newipla_error = ipl - ipl_ref_new;
  if ((activate != 0) && (enableBlock != 0)) {
    if (plateauDetect >= 1.0) {
      //  cq condition
      if ((plateauDetect - floor(plateauDetect / 2.0) * 2.0 == 0.0) &&
          (plateauDetectOld != plateauDetect)) {
        waitTime = time;
        enableChase = 1.0;
      }

      plateauDetectOld = plateauDetect;
      if ((time - waitTime >= delay_activate_ipl_cq) && (time - waitTimeChase >=
           delay_activate_ipl_chase) && (enableChase == 1.0) && (ipl_ref_new >=
           CQ_minIPL)) {
        guard1 = false;
        if (fabs(*newipla_error) >= threshold_ipl_drift_apart) {
          guard1 = true;
        } else {
          d0 = rt_roundd_snf(-vloop * (double)ipl_sign);
          if (d0 < 2.147483648E+9) {
            if (d0 >= -2.147483648E+9) {
              i0 = (int)d0;
            } else {
              i0 = MIN_int32_T;
            }
          } else if (d0 >= 2.147483648E+9) {
            i0 = MAX_int32_T;
          } else {
            i0 = 0;
          }

          if (i0 >= threshold_vloop) {
            guard1 = true;
          }
        }

        if (guard1) {
          if (plateauDetect - floor(plateauDetect / 2.0) * 2.0 == 0.0) {
            *enable = 1.0;
            newDelta = 0.0;
            if (fabs(*newipla_error) >= threshold_ipl_drift_apart) {
              newDelta = gamma_ipl_chase * *newipla_error;
            }

            // if vloop is due to mhd resistivity then
            d0 = rt_roundd_snf(-vloop * (double)ipl_sign);
            if (d0 < 2.147483648E+9) {
              if (d0 >= -2.147483648E+9) {
                i0 = (int)d0;
              } else {
                i0 = MIN_int32_T;
              }
            } else if (d0 >= 2.147483648E+9) {
              i0 = MAX_int32_T;
            } else {
              i0 = 0;
            }

            if ((i0 >= threshold_vloop) && (*newipla_error <=
                 threshold_ipl_vloop_activation)) {
              newDelta = -gamma_ipl_chase_vloop_mhd * ipl;
            }

            deltaIpl += newDelta;
            waitTimeChase = time;
          }
        }
      }
    }

    *deltaIpOut = deltaIpl;
  }
}

//
// Arguments    : void
// Return Type  : void
//
void ire_delta_ref_init()
{
  deltaIpl = 0.0;
  waitTime = 0.0;
  enableChase = 0.0;
  plateauDetectOld = 0.0;
  waitTimeChase = 0.0;
}

//
// File trailer for ire_delta_ref.cpp
//
// [EOF]
//

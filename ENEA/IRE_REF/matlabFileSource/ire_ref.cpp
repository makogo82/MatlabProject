//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ire_ref.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 15:51:22
//

// Include Files
#include "rt_nonfinite.h"
#include "ire_ref.h"

// Variable Definitions
static double enableCurrentDecay;
static double Iptpl;
static double tpl;
static double Ireftpl;
static double counter;

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
// % COMPUTE NEW REFERENCE FOR RUNAWAY ELECTRON CURRENT BEAM
//  INPUT:
//            time     =
//            activate = boolean var for eneabled the alghoritm
//                       0 : active
//                       1 : non active
//            ipl      = absolute value of the  plasma current
//            dipl     = derivative of abs(ipl)
//            ipl_ref  = prep_reference
//            plateauDetect =
//            deltaTdes     =
//            t_end         =
//            beta          =
//
//  OUTPUT:  y in R
//           y is the new reference for the ipl
// %
// Arguments    : double time
//                int activate
//                double ipl
//                double ipl_ref
//                int plateauDetect
//                double deltaTdes
//                double t_end
//                double beta
//                int enableMultiplePlateau
//                double percent_upgain_newiref
//                double threshold_time_newiref
//                double deltaIplChase
//                int enableUpdateRefUp
//                double *Iref
//                double *enable
//                double *statusChangeNewIref
// Return Type  : void
//
void ire_ref(double time, int activate, double ipl, double ipl_ref, int
             plateauDetect, double deltaTdes, double t_end, double beta, int
             enableMultiplePlateau, double percent_upgain_newiref, double
             threshold_time_newiref, double deltaIplChase, int enableUpdateRefUp,
             double *Iref, double *enable, double *statusChangeNewIref)
{
  int i0;
  double varargin_2;
  double b_deltaTdes;
  *Iref = ipl_ref;
  *enable = 0.0;
  *statusChangeNewIref = 0.0;
  if (activate != 0) {
    if ((plateauDetect >= 2) && (plateauDetect - ((plateauDetect >> 1) << 1) ==
         0)) {
      i0 = (int)rt_roundd_snf((double)plateauDetect / 2.0);
      if ((enableCurrentDecay == (double)i0 - 1.0) && ((plateauDetect == 2) ||
           (enableMultiplePlateau == 1))) {
        //  2 && 0, 4 && 1 , 6 && 2 ecc..
        // enableMultiplePlateau
        enableCurrentDecay++;

        // init var
        Iptpl = ipl;
        if (enableCurrentDecay == 1.0) {
          // only ones
          Ireftpl = ipl_ref;
        }

        tpl = time;
        counter = 0.0;
      }
    }

    if (enableCurrentDecay >= 1.0) {
      if (counter == 0.0) {
        if (enableCurrentDecay == 1.0) {
          Ireftpl = ipl_ref;
        }
      } else {
        // update program down
        // update program up
        if ((enableUpdateRefUp == 1) && (ipl - Ireftpl >= 0.0) && (t_end - time >=
             threshold_time_newiref)) {
          // update Iref_linear
          Iptpl = ipl + ipl * percent_upgain_newiref;
          *statusChangeNewIref = 1.0;
        }

        //  compute linear slop
        varargin_2 = t_end - tpl;

        //  compute exp decay
        if ((deltaTdes <= varargin_2) || rtIsNaN(varargin_2)) {
          b_deltaTdes = deltaTdes;
        } else {
          b_deltaTdes = varargin_2;
        }

        Ireftpl = beta * Ireftpl + (1.0 - beta) * ((Iptpl - Iptpl / b_deltaTdes *
          (time - tpl)) + deltaIplChase);
      }

      counter++;
      *Iref = Ireftpl;
      if (Ireftpl < 0.0) {
        //  disable
        enableCurrentDecay = -1.0;
      }
    }

    if (enableCurrentDecay == -1.0) {
      *Iref = 0.0;
    }

    *enable = enableCurrentDecay;
  }
}

//
// Arguments    : void
// Return Type  : void
//
void ire_ref_init()
{
  enableCurrentDecay = 0.0;
  Iptpl = 0.0;
  tpl = 0.0;
  Ireftpl = 0.0;
  counter = 0.0;
}

//
// File trailer for ire_ref.cpp
//
// [EOF]
//

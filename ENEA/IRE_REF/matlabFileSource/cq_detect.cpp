//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cq_detect.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 15:35:18
//

// Include Files
#include "rt_nonfinite.h"
#include "cq_detect.h"

// Variable Definitions
static double startDetection;

// Function Definitions

//
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
//                double ipl
//                double dipl
//                double minIPL
//                double dipl_threshold
//                double counterThreshold
// Return Type  : double
//
double cq_detect(int activate, double ipl, double dipl, double minIPL, double
                 dipl_threshold, double)
{
  double y;
  if (activate != 0) {
    // STEP1: FIND CQ
    if (ipl >= minIPL) {
      if ((dipl <= dipl_threshold) && (startDetection - floor(startDetection /
            2.0) * 2.0 == 0.0)) {
        startDetection++;
      }

      if ((dipl >= dipl_threshold) && (startDetection - floor(startDetection /
            2.0) * 2.0 == 1.0)) {
        startDetection++;
      }
    } else {
      if (startDetection > 0.0) {
        startDetection = -1.0;
      }
    }

    y = startDetection;
  } else {
    y = -2.0;
  }

  return y;
}

//
// Arguments    : void
// Return Type  : void
//
 

//
// File trailer for cq_detect.cpp
//
// [EOF]
//

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cq_detect.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 15:35:18
//
#ifndef CQ_DETECT_H
#define CQ_DETECT_H

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "cq_detect_types.h"

// Function Declarations
extern double cq_detect(int activate, double ipl, double dipl, double minIPL,
  double dipl_threshold, double counterThreshold);
extern void cq_detect_init();

#endif

//
// File trailer for cq_detect.h
//
// [EOF]
//

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ire_delta_ref.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 16:09:20
//
#ifndef IRE_DELTA_REF_H
#define IRE_DELTA_REF_H

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "ire_delta_ref_types.h"

// Function Declarations
extern void ire_delta_ref(double time, int activate, double vloop, double ipl,
  int ipl_sign, double ipl_ref_new, double plateauDetect, double
  threshold_ipl_drift_apart, double gamma_ipl_chase, int
  delay_activate_ipl_chase, double delay_activate_ipl_cq, double threshold_vloop,
  double threshold_ipl_vloop_activation, double gamma_ipl_chase_vloop_mhd, int
  enableBlock, double CQ_minIPL, double *deltaIpOut, double *newipla_error,
  double *enable);
extern void ire_delta_ref_init();

#endif

//
// File trailer for ire_delta_ref.h
//
// [EOF]
//

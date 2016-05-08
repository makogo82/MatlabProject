//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: square_wave.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 16:17:07
//
#ifndef SQUARE_WAVE_H
#define SQUARE_WAVE_H

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "square_wave_types.h"

// Function Declarations
extern void square_wave(double time, int activate, double ipl_ref, double
  ipl_ref_new, double CQ_minIPL, double plateauDetect, double pwm_amp, double
  pwm_period, double delay_activate_pwm, double duty_cycle_pwm, int enableBlock,
  double *deltaIpOut, double *enable);
extern void square_wave_init();

#endif

//
// File trailer for square_wave.h
//
// [EOF]
//

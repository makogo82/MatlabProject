//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ire_ref.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Apr-2016 15:51:22
//
#ifndef IRE_REF_H
#define IRE_REF_H

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "ire_ref_types.h"

// Function Declarations
extern void ire_ref(double time, int activate, double ipl, double ipl_ref, int
                    plateauDetect, double deltaTdes, double t_end, double beta,
                    int enableMultiplePlateau, double percent_upgain_newiref,
                    double threshold_time_newiref, double deltaIplChase, int
                    enableUpdateRefUp, double *Iref, double *enable, double
                    *statusChangeNewIref);
extern void ire_ref_init();

#endif

//
// File trailer for ire_ref.h
//
// [EOF]
//

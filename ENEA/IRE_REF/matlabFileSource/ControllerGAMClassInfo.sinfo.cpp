#define protected public
#define private public
#include "ControllerGAMClassInfo.h"
#include "ObjectRegistryItem.h"
#include "ClassStructure.h"
#include "ObjectMacros.h"
static ClassStructureEntry ControllerGAMOutputStructure_IfControl_CSE_EL("float","",0,0,0,0,0 ,"IfControl",msizeof(ControllerGAMOutputStructure,IfControl),indexof(ControllerGAMOutputStructure,IfControl));
static ClassStructureEntry ControllerGAMOutputStructure_IvControl_CSE_EL("float","",0,0,0,0,0 ,"IvControl",msizeof(ControllerGAMOutputStructure,IvControl),indexof(ControllerGAMOutputStructure,IvControl));
static ClassStructureEntry ControllerGAMOutputStructure_IhControl_CSE_EL("float","",0,0,0,0,0 ,"IhControl",msizeof(ControllerGAMOutputStructure,IhControl),indexof(ControllerGAMOutputStructure,IhControl));
static ClassStructureEntry ControllerGAMOutputStructure_ItControl_CSE_EL("float","",0,0,0,0,0 ,"ItControl",msizeof(ControllerGAMOutputStructure,ItControl),indexof(ControllerGAMOutputStructure,ItControl));
static ClassStructureEntry ControllerGAMOutputStructure_IPlasmaOn_CSE_EL("float","",0,0,0,0,0 ,"IPlasmaOn",msizeof(ControllerGAMOutputStructure,IPlasmaOn),indexof(ControllerGAMOutputStructure,IPlasmaOn));
static ClassStructureEntry ControllerGAMOutputStructure_awv1_CSE_EL("float","",0,0,0,0,0 ,"awv1",msizeof(ControllerGAMOutputStructure,awv1),indexof(ControllerGAMOutputStructure,awv1));
static ClassStructureEntry ControllerGAMOutputStructure_awv2_CSE_EL("float","",0,0,0,0,0 ,"awv2",msizeof(ControllerGAMOutputStructure,awv2),indexof(ControllerGAMOutputStructure,awv2));
static ClassStructureEntry ControllerGAMOutputStructure_allocDeltaIf_CSE_EL("float","",0,0,0,0,0 ,"allocDeltaIf",msizeof(ControllerGAMOutputStructure,allocDeltaIf),indexof(ControllerGAMOutputStructure,allocDeltaIf));
static ClassStructureEntry ControllerGAMOutputStructure_allocDeltaIv_CSE_EL("float","",0,0,0,0,0 ,"allocDeltaIv",msizeof(ControllerGAMOutputStructure,allocDeltaIv),indexof(ControllerGAMOutputStructure,allocDeltaIv));
static ClassStructureEntry ControllerGAMOutputStructure_allocIFRef_CSE_EL("float","",0,0,0,0,0 ,"allocIFRef",msizeof(ControllerGAMOutputStructure,allocIFRef),indexof(ControllerGAMOutputStructure,allocIFRef));
static ClassStructureEntry ControllerGAMOutputStructure_allocVerticalFieldBalancer_CSE_EL("float","",0,0,0,0,0 ,"allocVerticalFieldBalancer",msizeof(ControllerGAMOutputStructure,allocVerticalFieldBalancer),indexof(ControllerGAMOutputStructure,allocVerticalFieldBalancer));
static ClassStructureEntry ControllerGAMOutputStructure_f1_CSE_EL("float","",0,0,0,0,0 ,"f1",msizeof(ControllerGAMOutputStructure,f1),indexof(ControllerGAMOutputStructure,f1));
static ClassStructureEntry ControllerGAMOutputStructure_f2_CSE_EL("float","",0,0,0,0,0 ,"f2",msizeof(ControllerGAMOutputStructure,f2),indexof(ControllerGAMOutputStructure,f2));
static ClassStructureEntry ControllerGAMOutputStructure_f3_CSE_EL("float","",0,0,0,0,0 ,"f3",msizeof(ControllerGAMOutputStructure,f3),indexof(ControllerGAMOutputStructure,f3));
static ClassStructureEntry ControllerGAMOutputStructure_f4_CSE_EL("float","",0,0,0,0,0 ,"f4",msizeof(ControllerGAMOutputStructure,f4),indexof(ControllerGAMOutputStructure,f4));
static ClassStructureEntry ControllerGAMOutputStructure_f5_CSE_EL("float","",0,0,0,0,0 ,"f5",msizeof(ControllerGAMOutputStructure,f5),indexof(ControllerGAMOutputStructure,f5));
static ClassStructureEntry ControllerGAMOutputStructure_f6_CSE_EL("float","",0,0,0,0,0 ,"f6",msizeof(ControllerGAMOutputStructure,f6),indexof(ControllerGAMOutputStructure,f6));
static ClassStructureEntry ControllerGAMOutputStructure_f7_CSE_EL("float","",0,0,0,0,0 ,"f7",msizeof(ControllerGAMOutputStructure,f7),indexof(ControllerGAMOutputStructure,f7));
static ClassStructureEntry ControllerGAMOutputStructure_kfabort_CSE_EL("float","",0,0,0,0,0 ,"kfabort",msizeof(ControllerGAMOutputStructure,kfabort),indexof(ControllerGAMOutputStructure,kfabort));
static ClassStructureEntry ControllerGAMOutputStructure_kvtabort_CSE_EL("float","",0,0,0,0,0 ,"kvtabort",msizeof(ControllerGAMOutputStructure,kvtabort),indexof(ControllerGAMOutputStructure,kvtabort));
static ClassStructureEntry ControllerGAMOutputStructure_fl_plasma_CSE_EL("int","",0,0,0,0,0 ,"fl_plasma",msizeof(ControllerGAMOutputStructure,fl_plasma),indexof(ControllerGAMOutputStructure,fl_plasma));
static ClassStructureEntry ControllerGAMOutputStructure_kpsf_CSE_EL("float","",0,0,0,0,0 ,"kpsf",msizeof(ControllerGAMOutputStructure,kpsf),indexof(ControllerGAMOutputStructure,kpsf));
static ClassStructureEntry ControllerGAMOutputStructure_kisf_CSE_EL("float","",0,0,0,0,0 ,"kisf",msizeof(ControllerGAMOutputStructure,kisf),indexof(ControllerGAMOutputStructure,kisf));
static ClassStructureEntry ControllerGAMOutputStructure_kdsf_CSE_EL("float","",0,0,0,0,0 ,"kdsf",msizeof(ControllerGAMOutputStructure,kdsf),indexof(ControllerGAMOutputStructure,kdsf));
static ClassStructureEntry ControllerGAMOutputStructure_pidF_CSE_EL("float","",0,0,0,0,0 ,"pidF",msizeof(ControllerGAMOutputStructure,pidF),indexof(ControllerGAMOutputStructure,pidF));
static ClassStructureEntry ControllerGAMOutputStructure_pidF2_CSE_EL("float","",0,0,0,0,0 ,"pidF2",msizeof(ControllerGAMOutputStructure,pidF2),indexof(ControllerGAMOutputStructure,pidF2));
static ClassStructureEntry ControllerGAMOutputStructure_exp_CSE_EL("float","",0,0,0,0,0 ,"exp",msizeof(ControllerGAMOutputStructure,exp),indexof(ControllerGAMOutputStructure,exp));
static ClassStructureEntry ControllerGAMOutputStructure_e_hat_CSE_EL("float","",0,0,0,0,0 ,"e_hat",msizeof(ControllerGAMOutputStructure,e_hat),indexof(ControllerGAMOutputStructure,e_hat));
static ClassStructureEntry ControllerGAMOutputStructure_xi_CSE_EL("float","",0,0,0,0,0 ,"xi",msizeof(ControllerGAMOutputStructure,xi),indexof(ControllerGAMOutputStructure,xi));
static ClassStructureEntry ControllerGAMOutputStructure_et_CSE_EL("float","",0,0,0,0,0 ,"et",msizeof(ControllerGAMOutputStructure,et),indexof(ControllerGAMOutputStructure,et));
static ClassStructureEntry ControllerGAMOutputStructure_altnew_CSE_EL("float","",0,0,0,0,0 ,"altnew",msizeof(ControllerGAMOutputStructure,altnew),indexof(ControllerGAMOutputStructure,altnew));
static ClassStructureEntry ControllerGAMOutputStructure_xir_CSE_EL("float","",0,0,0,0,0 ,"xir",msizeof(ControllerGAMOutputStructure,xir),indexof(ControllerGAMOutputStructure,xir));
static ClassStructureEntry ControllerGAMOutputStructure_altnewr_CSE_EL("float","",0,0,0,0,0 ,"altnewr",msizeof(ControllerGAMOutputStructure,altnewr),indexof(ControllerGAMOutputStructure,altnewr));
static ClassStructureEntry ControllerGAMOutputStructure_ZFMRamp_CSE_EL("float","",0,0,0,0,0 ,"ZFMRamp",msizeof(ControllerGAMOutputStructure,ZFMRamp),indexof(ControllerGAMOutputStructure,ZFMRamp));
static ClassStructureEntry ControllerGAMOutputStructure_ZFM_I_C_CSE_EL("float","",0,0,0,0,0 ,"ZFM_I_C",msizeof(ControllerGAMOutputStructure,ZFM_I_C),indexof(ControllerGAMOutputStructure,ZFM_I_C));
static ClassStructureEntry ControllerGAMOutputStructure_alfnew_CSE_EL("float","",0,0,0,0,0 ,"alfnew",msizeof(ControllerGAMOutputStructure,alfnew),indexof(ControllerGAMOutputStructure,alfnew));
static ClassStructureEntry ControllerGAMOutputStructure_dIp_CSE_EL("float","",0,0,0,0,0 ,"dIp",msizeof(ControllerGAMOutputStructure,dIp),indexof(ControllerGAMOutputStructure,dIp));
static ClassStructureEntry ControllerGAMOutputStructure_ipl_sign_CSE_EL("float","",0,0,0,0,0 ,"ipl_sign",msizeof(ControllerGAMOutputStructure,ipl_sign),indexof(ControllerGAMOutputStructure,ipl_sign));
static ClassStructureEntry ControllerGAMOutputStructure_cq_detect_output_CSE_EL("float","",0,0,0,0,0 ,"cq_detect_output",msizeof(ControllerGAMOutputStructure,cq_detect_output),indexof(ControllerGAMOutputStructure,cq_detect_output));
static ClassStructureEntry ControllerGAMOutputStructure_iref_output_CSE_EL("float","",0,0,0,0,0 ,"iref_output",msizeof(ControllerGAMOutputStructure,iref_output),indexof(ControllerGAMOutputStructure,iref_output));
static ClassStructureEntry ControllerGAMOutputStructure_deltaIpOut_CSE_EL("float","",0,0,0,0,0 ,"deltaIpOut",msizeof(ControllerGAMOutputStructure,deltaIpOut),indexof(ControllerGAMOutputStructure,deltaIpOut));
static ClassStructureEntry ControllerGAMOutputStructure_deltaIpOut_pwm_CSE_EL("float","",0,0,0,0,0 ,"deltaIpOut_pwm",msizeof(ControllerGAMOutputStructure,deltaIpOut_pwm),indexof(ControllerGAMOutputStructure,deltaIpOut_pwm));
static ClassStructureEntry ControllerGAMOutputStructure_enable_pwm_CSE_EL("int","",0,0,0,0,0 ,"enable_pwm",msizeof(ControllerGAMOutputStructure,enable_pwm),indexof(ControllerGAMOutputStructure,enable_pwm));
static ClassStructureEntry ControllerGAMOutputStructure_enable_delta_ref_CSE_EL("int","",0,0,0,0,0 ,"enable_delta_ref",msizeof(ControllerGAMOutputStructure,enable_delta_ref),indexof(ControllerGAMOutputStructure,enable_delta_ref));
static ClassStructureEntry ControllerGAMOutputStructure_enable_current_tracking_CSE_EL("int","",0,0,0,0,0 ,"enable_current_tracking",msizeof(ControllerGAMOutputStructure,enable_current_tracking),indexof(ControllerGAMOutputStructure,enable_current_tracking));
static ClassStructureEntry ControllerGAMOutputStructure_ipl_error_CSE_EL("float","",0,0,0,0,0 ,"ipl_error",msizeof(ControllerGAMOutputStructure,ipl_error),indexof(ControllerGAMOutputStructure,ipl_error));
static ClassStructureEntry ControllerGAMOutputStructure_ipl_ref_old_CSE_EL("float","",0,0,0,0,0 ,"ipl_ref_old",msizeof(ControllerGAMOutputStructure,ipl_ref_old),indexof(ControllerGAMOutputStructure,ipl_ref_old));
static ClassStructureEntry ControllerGAMOutputStructure_time_CSE_EL("float","",0,0,0,0,0 ,"time",msizeof(ControllerGAMOutputStructure,time),indexof(ControllerGAMOutputStructure,time));
static ClassStructureEntry ControllerGAMOutputStructure_iref_linear_CSE_EL("float","",0,0,0,0,0 ,"iref_linear",msizeof(ControllerGAMOutputStructure,iref_linear),indexof(ControllerGAMOutputStructure,iref_linear));
static ClassStructureEntry * ControllerGAMOutputStructure__CSE__[] = {
    &ControllerGAMOutputStructure_IfControl_CSE_EL,
    &ControllerGAMOutputStructure_IvControl_CSE_EL,
    &ControllerGAMOutputStructure_IhControl_CSE_EL,
    &ControllerGAMOutputStructure_ItControl_CSE_EL,
    &ControllerGAMOutputStructure_IPlasmaOn_CSE_EL,
    &ControllerGAMOutputStructure_awv1_CSE_EL,
    &ControllerGAMOutputStructure_awv2_CSE_EL,
    &ControllerGAMOutputStructure_allocDeltaIf_CSE_EL,
    &ControllerGAMOutputStructure_allocDeltaIv_CSE_EL,
    &ControllerGAMOutputStructure_allocIFRef_CSE_EL,
    &ControllerGAMOutputStructure_allocVerticalFieldBalancer_CSE_EL,
    &ControllerGAMOutputStructure_f1_CSE_EL,
    &ControllerGAMOutputStructure_f2_CSE_EL,
    &ControllerGAMOutputStructure_f3_CSE_EL,
    &ControllerGAMOutputStructure_f4_CSE_EL,
    &ControllerGAMOutputStructure_f5_CSE_EL,
    &ControllerGAMOutputStructure_f6_CSE_EL,
    &ControllerGAMOutputStructure_f7_CSE_EL,
    &ControllerGAMOutputStructure_kfabort_CSE_EL,
    &ControllerGAMOutputStructure_kvtabort_CSE_EL,
    &ControllerGAMOutputStructure_fl_plasma_CSE_EL,
    &ControllerGAMOutputStructure_kpsf_CSE_EL,
    &ControllerGAMOutputStructure_kisf_CSE_EL,
    &ControllerGAMOutputStructure_kdsf_CSE_EL,
    &ControllerGAMOutputStructure_pidF_CSE_EL,
    &ControllerGAMOutputStructure_pidF2_CSE_EL,
    &ControllerGAMOutputStructure_exp_CSE_EL,
    &ControllerGAMOutputStructure_e_hat_CSE_EL,
    &ControllerGAMOutputStructure_xi_CSE_EL,
    &ControllerGAMOutputStructure_et_CSE_EL,
    &ControllerGAMOutputStructure_altnew_CSE_EL,
    &ControllerGAMOutputStructure_xir_CSE_EL,
    &ControllerGAMOutputStructure_altnewr_CSE_EL,
    &ControllerGAMOutputStructure_ZFMRamp_CSE_EL,
    &ControllerGAMOutputStructure_ZFM_I_C_CSE_EL,
    &ControllerGAMOutputStructure_alfnew_CSE_EL,
    &ControllerGAMOutputStructure_dIp_CSE_EL,
    &ControllerGAMOutputStructure_ipl_sign_CSE_EL,
    &ControllerGAMOutputStructure_cq_detect_output_CSE_EL,
    &ControllerGAMOutputStructure_iref_output_CSE_EL,
    &ControllerGAMOutputStructure_deltaIpOut_CSE_EL,
    &ControllerGAMOutputStructure_deltaIpOut_pwm_CSE_EL,
    &ControllerGAMOutputStructure_enable_pwm_CSE_EL,
    &ControllerGAMOutputStructure_enable_delta_ref_CSE_EL,
    &ControllerGAMOutputStructure_enable_current_tracking_CSE_EL,
    &ControllerGAMOutputStructure_ipl_error_CSE_EL,
    &ControllerGAMOutputStructure_ipl_ref_old_CSE_EL,
    &ControllerGAMOutputStructure_time_CSE_EL,
    &ControllerGAMOutputStructure_iref_linear_CSE_EL,
    NULL
};
ClassStructure ControllerGAMOutputStructure__CS__("ControllerGAMOutputStructure",sizeof(ControllerGAMOutputStructure),0 ,ControllerGAMOutputStructure__CSE__);
STRUCTREGISTER("ControllerGAMOutputStructure",ControllerGAMOutputStructure__CS__)
static ClassStructureEntry ControllerGAMInputStructure_usecTime_CSE_EL("int","",0,0,0,0,0 ,"usecTime",msizeof(ControllerGAMInputStructure,usecTime),indexof(ControllerGAMInputStructure,usecTime));
static ClassStructureEntry ControllerGAMInputStructure_fbcor_ALT_prep_CSE_EL("float","",0,0,0,0,0 ,"fbcor_ALT_prep",msizeof(ControllerGAMInputStructure,fbcor_ALT_prep),indexof(ControllerGAMInputStructure,fbcor_ALT_prep));
static ClassStructureEntry ControllerGAMInputStructure_fbcor_ALF_prep_CSE_EL("float","",0,0,0,0,0 ,"fbcor_ALF_prep",msizeof(ControllerGAMInputStructure,fbcor_ALF_prep),indexof(ControllerGAMInputStructure,fbcor_ALF_prep));
static ClassStructureEntry ControllerGAMInputStructure_fbcor_ALV_prep_CSE_EL("float","",0,0,0,0,0 ,"fbcor_ALV_prep",msizeof(ControllerGAMInputStructure,fbcor_ALV_prep),indexof(ControllerGAMInputStructure,fbcor_ALV_prep));
static ClassStructureEntry ControllerGAMInputStructure_fbcor_ALH_prep_CSE_EL("float","",0,0,0,0,0 ,"fbcor_ALH_prep",msizeof(ControllerGAMInputStructure,fbcor_ALH_prep),indexof(ControllerGAMInputStructure,fbcor_ALH_prep));
static ClassStructureEntry ControllerGAMInputStructure_ipl_CSE_EL("float","",0,0,0,0,0 ,"ipl",msizeof(ControllerGAMInputStructure,ipl),indexof(ControllerGAMInputStructure,ipl));
static ClassStructureEntry ControllerGAMInputStructure_ipRif_CSE_EL("float","",0,0,0,0,0 ,"ipRif",msizeof(ControllerGAMInputStructure,ipRif),indexof(ControllerGAMInputStructure,ipRif));
static ClassStructureEntry ControllerGAMInputStructure_ipRif1_CSE_EL("float","",0,0,0,0,0 ,"ipRif1",msizeof(ControllerGAMInputStructure,ipRif1),indexof(ControllerGAMInputStructure,ipRif1));
static ClassStructureEntry ControllerGAMInputStructure_dep_CSE_EL("float","",0,0,0,0,0 ,"dep",msizeof(ControllerGAMInputStructure,dep),indexof(ControllerGAMInputStructure,dep));
static ClassStructureEntry ControllerGAMInputStructure_dez_CSE_EL("float","",0,0,0,0,0 ,"dez",msizeof(ControllerGAMInputStructure,dez),indexof(ControllerGAMInputStructure,dez));
static ClassStructureEntry ControllerGAMInputStructure_fbTenVpRif_CSE_EL("float","",0,0,0,0,0 ,"fbTenVpRif",msizeof(ControllerGAMInputStructure,fbTenVpRif),indexof(ControllerGAMInputStructure,fbTenVpRif));
static ClassStructureEntry ControllerGAMInputStructure_fbadmBpRif_CSE_EL("float","",0,0,0,0,0 ,"fbadmBpRif",msizeof(ControllerGAMInputStructure,fbadmBpRif),indexof(ControllerGAMInputStructure,fbadmBpRif));
static ClassStructureEntry ControllerGAMInputStructure_fbindLipRif_CSE_EL("float","",0,0,0,0,0 ,"fbindLipRif",msizeof(ControllerGAMInputStructure,fbindLipRif),indexof(ControllerGAMInputStructure,fbindLipRif));
static ClassStructureEntry ControllerGAMInputStructure_badGasMeas_CSE_EL("int","",0,0,0,0,0 ,"badGasMeas",msizeof(ControllerGAMInputStructure,badGasMeas),indexof(ControllerGAMInputStructure,badGasMeas));
static ClassStructureEntry ControllerGAMInputStructure_runaway_CSE_EL("int","",0,0,0,0,0 ,"runaway",msizeof(ControllerGAMInputStructure,runaway),indexof(ControllerGAMInputStructure,runaway));
static ClassStructureEntry ControllerGAMInputStructure_disruption_CSE_EL("int","",0,0,0,0,0 ,"disruption",msizeof(ControllerGAMInputStructure,disruption),indexof(ControllerGAMInputStructure,disruption));
static ClassStructureEntry ControllerGAMInputStructure_runawayPlateau_CSE_EL("int","",0,0,0,0,0 ,"runawayPlateau",msizeof(ControllerGAMInputStructure,runawayPlateau),indexof(ControllerGAMInputStructure,runawayPlateau));
static ClassStructureEntry ControllerGAMInputStructure_elong_CSE_EL("float","",0,0,0,0,0 ,"elong",msizeof(ControllerGAMInputStructure,elong),indexof(ControllerGAMInputStructure,elong));
static ClassStructureEntry ControllerGAMInputStructure_slowDeIP_CSE_EL("float","",0,0,0,0,0 ,"slowDeIP",msizeof(ControllerGAMInputStructure,slowDeIP),indexof(ControllerGAMInputStructure,slowDeIP));
static ClassStructureEntry ControllerGAMInputStructure_f8_CSE_EL("float","",0,0,0,0,0 ,"f8",msizeof(ControllerGAMInputStructure,f8),indexof(ControllerGAMInputStructure,f8));
static ClassStructureEntry ControllerGAMInputStructure_ZTM_I_CSE_EL("float","",0,0,0,0,0 ,"ZTM_I",msizeof(ControllerGAMInputStructure,ZTM_I),indexof(ControllerGAMInputStructure,ZTM_I));
static ClassStructureEntry ControllerGAMInputStructure_e_hat_CSE_EL("float","",0,0,0,0,0 ,"e_hat",msizeof(ControllerGAMInputStructure,e_hat),indexof(ControllerGAMInputStructure,e_hat));
static ClassStructureEntry ControllerGAMInputStructure_VLOOP_CSE_EL("float","",0,0,0,0,0 ,"VLOOP",msizeof(ControllerGAMInputStructure,VLOOP),indexof(ControllerGAMInputStructure,VLOOP));
static ClassStructureEntry ControllerGAMInputStructure_deVLOOP_CSE_EL("float","",0,0,0,0,0 ,"deVLOOP",msizeof(ControllerGAMInputStructure,deVLOOP),indexof(ControllerGAMInputStructure,deVLOOP));
static ClassStructureEntry ControllerGAMInputStructure_xir_CSE_EL("float","",0,0,0,0,0 ,"xir",msizeof(ControllerGAMInputStructure,xir),indexof(ControllerGAMInputStructure,xir));
static ClassStructureEntry ControllerGAMInputStructure_ZFM_I_CSE_EL("float","",0,0,0,0,0 ,"ZFM_I",msizeof(ControllerGAMInputStructure,ZFM_I),indexof(ControllerGAMInputStructure,ZFM_I));
static ClassStructureEntry * ControllerGAMInputStructure__CSE__[] = {
    &ControllerGAMInputStructure_usecTime_CSE_EL,
    &ControllerGAMInputStructure_fbcor_ALT_prep_CSE_EL,
    &ControllerGAMInputStructure_fbcor_ALF_prep_CSE_EL,
    &ControllerGAMInputStructure_fbcor_ALV_prep_CSE_EL,
    &ControllerGAMInputStructure_fbcor_ALH_prep_CSE_EL,
    &ControllerGAMInputStructure_ipl_CSE_EL,
    &ControllerGAMInputStructure_ipRif_CSE_EL,
    &ControllerGAMInputStructure_ipRif1_CSE_EL,
    &ControllerGAMInputStructure_dep_CSE_EL,
    &ControllerGAMInputStructure_dez_CSE_EL,
    &ControllerGAMInputStructure_fbTenVpRif_CSE_EL,
    &ControllerGAMInputStructure_fbadmBpRif_CSE_EL,
    &ControllerGAMInputStructure_fbindLipRif_CSE_EL,
    &ControllerGAMInputStructure_badGasMeas_CSE_EL,
    &ControllerGAMInputStructure_runaway_CSE_EL,
    &ControllerGAMInputStructure_disruption_CSE_EL,
    &ControllerGAMInputStructure_runawayPlateau_CSE_EL,
    &ControllerGAMInputStructure_elong_CSE_EL,
    &ControllerGAMInputStructure_slowDeIP_CSE_EL,
    &ControllerGAMInputStructure_f8_CSE_EL,
    &ControllerGAMInputStructure_ZTM_I_CSE_EL,
    &ControllerGAMInputStructure_e_hat_CSE_EL,
    &ControllerGAMInputStructure_VLOOP_CSE_EL,
    &ControllerGAMInputStructure_deVLOOP_CSE_EL,
    &ControllerGAMInputStructure_xir_CSE_EL,
    &ControllerGAMInputStructure_ZFM_I_CSE_EL,
    NULL
};
ClassStructure ControllerGAMInputStructure__CS__("ControllerGAMInputStructure",sizeof(ControllerGAMInputStructure),0 ,ControllerGAMInputStructure__CSE__);
STRUCTREGISTER("ControllerGAMInputStructure",ControllerGAMInputStructure__CS__)
static ClassStructureEntry ControllerGAMClassInfo_input_CSE_EL("ControllerGAMInputStructure","",0,0,0,0,0 ,"input",msizeof(ControllerGAMClassInfo,input),indexof(ControllerGAMClassInfo,input));
static ClassStructureEntry ControllerGAMClassInfo_output_CSE_EL("ControllerGAMOutputStructure","",0,0,0,0,0 ,"output",msizeof(ControllerGAMClassInfo,output),indexof(ControllerGAMClassInfo,output));
static ClassStructureEntry * ControllerGAMClassInfo__CSE__[] = {
    &ControllerGAMClassInfo_input_CSE_EL,
    &ControllerGAMClassInfo_output_CSE_EL,
    NULL
};
ClassStructure ControllerGAMClassInfo__CS__("ControllerGAMClassInfo",sizeof(ControllerGAMClassInfo),0 ,ControllerGAMClassInfo__CSE__);
STRUCTREGISTER("ControllerGAMClassInfo",ControllerGAMClassInfo__CS__)
ClassStructure * ControllerGAMClassInfo_sinfo[] = {
    &ControllerGAMOutputStructure__CS__,
    &ControllerGAMInputStructure__CS__,
    &ControllerGAMClassInfo__CS__,
    NULL
};
/* All user tunable parameters, including functional weights,
   should be listed here.

   the PARAM macro to generates both and enums and strings from the
   parameter names. Add new parameters at the end of this list, but
   keep XC_NR_PARAMS last  */

  PARAM(XC_SLATERX),
  PARAM(XC_VWN5C),
  PARAM(XC_BECKEX),
  PARAM(XC_BECKECORRX),
  PARAM(XC_BECKESRX),
  PARAM(XC_LB94),
  PARAM(XC_LYPC),
  PARAM(XC_PBEX),
  PARAM(XC_REVPBEX),
  PARAM(XC_PBEC),
  PARAM(XC_LDAERFX),
  PARAM(XC_LDAERFC),
  PARAM(XC_RANGESEP_MU),
  PARAM(XC_KTX),
  PARAM(XC_TFK),
  PARAM(XC_PW91X),
  PARAM(XC_PW91K),
  PARAM(XC_PW92C),
  PARAM(XC_M05X),
  PARAM(XC_M05X2X),
  PARAM(XC_M06X),
  PARAM(XC_M06X2X),
  PARAM(XC_M06LX),
  PARAM(XC_M06HFX),
  // In development
  PARAM(XC_M05X2C),
  PARAM(XC_M05C),
  PARAM(XC_M06C),
  PARAM(XC_M06X2C),
  // This has to be last:
  PARAM(XC_NR_PARAMS) 

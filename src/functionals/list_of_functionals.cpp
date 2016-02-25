#include "functional.h"

/* See also list_of_parameters.h for list of functional
   parameters. */

#define SETUP(F) void F(functional &); xc_run_functional_setup(F);

void xcint_setup_functionals()
{
  SETUP(setup_slaterx);
  SETUP(setup_vwn5);
  SETUP(setup_beckex);
  SETUP(setup_beckexcorr);
  SETUP(setup_beckexsr);
  SETUP(setup_optx);
  //  SETUP(setup_lb94);
  SETUP(setup_lypc);

  SETUP(setup_pbex);
  SETUP(setup_revpbex);
  SETUP(setup_rpbex);
  SETUP(setup_pbec);
  SETUP(setup_spbec);
  SETUP(setup_vwn_pbec);
#ifndef XCFUN_NO_ERF
  SETUP(setup_ldaerfx);
  SETUP(setup_ldaerfc);
  SETUP(setup_ldaerfc_jt);
#endif
  SETUP(setup_ktx);
  SETUP(setup_tfk);
  SETUP(setup_pw91x);
  SETUP(setup_pw91k);
  SETUP(setup_pw92c);
  SETUP(setup_m05x);
  SETUP(setup_m05x2x);
  SETUP(setup_m06x);
  SETUP(setup_m06x2x);
  SETUP(setup_m06lx);
  SETUP(setup_m06hfx);
#ifdef XCFUN_IN_DEVELOPMENT
  SETUP(setup_m05x2c);
  SETUP(setup_m05c);
  SETUP(setup_m06c);
  SETUP(setup_m06hfc);
  SETUP(setup_m06lc);
  SETUP(setup_m06x2c);
  SETUP(setup_tpssc);
  SETUP(setup_tpssx);
  SETUP(setup_revtpssc);
  SETUP(setup_revtpssx);
#endif
  SETUP(setup_brx);
}

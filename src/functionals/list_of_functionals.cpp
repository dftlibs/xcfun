#include "functional.h"

#define SETUP(F) void F(functional &); xc_run_functional_setup(F);

void xcint_setup_functionals()
{
  SETUP(setup_slaterx);
  SETUP(setup_vwn5);
  SETUP(setup_beckex);
  SETUP(setup_beckexcorr);
  SETUP(setup_beckexsr);
  SETUP(setup_lb94);
  SETUP(setup_lypc);
  SETUP(setup_pbex);
  SETUP(setup_pbec);
  SETUP(setup_ldaerfx);
  SETUP(setup_ldaerfc);
  SETUP(setup_ktx);
  SETUP(setup_tfk);
  SETUP(setup_pw91x);
  SETUP(setup_pw91k);
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
  SETUP(setup_m06x2c);
#endif
}

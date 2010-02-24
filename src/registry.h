#ifndef REGISTRY_H
#define REGISTRY_H
#include <ostream>
#include "densvars.h"
#include "xc_fun.h"

template<enum xc_funparams functional_id>
struct functional
{
  static char *get_name(void);
  static char *get_reference(void);
  template<class num>
  static num energy(const densvars<num> &d);
  static int test(void);
};

template<>
struct functional<XC_SLATER_EXCHANGE>
{
  static const char *get_name(void) { return "slaterx"; }
  static const char *get_reference(void) { return "blah blah"; }
  template<class num>
  static num energy(const densvars<num> &d) { return -12; }
  static int test(void) { return 0; }
};




#endif

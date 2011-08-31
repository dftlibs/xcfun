#ifndef SLATER_H
#define SLATER_H

#include "constants.hpp"

template<class num>
static num slaterx(const densvars<num> &d) 
{ 
  return (-xc_constants::c_slater)*(d.a_43 + d.b_43);
}

#endif

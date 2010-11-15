#include <cstdlib>
#include "functional.h"
#include "xcint.h"

template<class T>
double densvars<T>::get_param(enum xc_parameter p) const
{
  return parent->settings[p];
}

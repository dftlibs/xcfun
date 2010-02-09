// This file is for constants that are used in more than one
// functional. It will be included exactly once, by xc_fun.cpp.
// Don't include it in your functional .h file. 

#ifndef M_PI // M_PI is not standard for some reason
#define M_PI 3.14159265358979323846
#endif

namespace xc_constants
{
  static const double CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
};


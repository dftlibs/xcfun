#ifndef SLATER_H
#define SLATER_H

template<class num>
static num slaterx(const densvars<num> &d) 
{ 
  return -3.0/4.0*pow(6/M_PI, 1.0/3.0)*
    (pow(d.a,4.0/3.0) + pow(d.b,4.0/3.0));
}

#endif

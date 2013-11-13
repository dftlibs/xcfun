#include "functional.hpp"
#include "constants.hpp"

template<class num>
static num b97c(const densvars<num> &d) 
{
  return 0;
}

template<class num>
static num b97x(const densvars<num> &d) 
{
  return 0;
}

FUNCTIONAL(XC_B97C) = {
  "B97 correlation",
  "reference goes here",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(b97c) };


FUNCTIONAL(XC_B97X) = {
  "B97 exchange",
  "reference goes here",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(b97x) };

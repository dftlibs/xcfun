#include "functional.hpp"
#include "constants.hpp"
#include "b97x.hpp"
#include "b97c.hpp"


template<class num>
static num b97x_en(const densvars<num> &d)
{
    return b97x::energy_b97x_ab(b97x::Gamma,b97x::c_b97,d.a_43,d.gaa)
    + b97x::energy_b97x_ab(b97x::Gamma,b97x::c_b97, d.b_43,d.gbb);
}


template<class num>
static num b97c_en(const densvars<num> &d)
{
    num e_LSDA_a, e_LSDA_b, tmp;
    tmp = b97c::energy_b97c_par(b97c::Gamma_par, b97c::c_b97[1] ,
                                d.a,d.a_43, d.gaa, e_LSDA_a) +
    b97c::energy_b97c_par(b97c::Gamma_par, b97c::c_b97[1] , d.b,  d.b_43, d.gbb, e_LSDA_b);
    
    return tmp + b97c::energy_b97c_antipar(b97c::Gamma_antipar,
                                           b97c::c_b97[0], d, e_LSDA_a, e_LSDA_b);

}


FUNCTIONAL(XC_B97X) = {
  "B97 exchange",
  "Density-functional thermochemistry. V.\n"
  "Systematic optimization of exchange-correlation functionals\n"
  "A.D.Becke; J.Chem.Phys.; 107, 8554, (1997)\n"
  "Implemented by Sarah Reimann ",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(b97x_en) };



FUNCTIONAL(XC_B97C) = {
  "B97 correlation",
    "Density-functional thermochemistry\n"
    "Systematic optimization of exchange-correlation functionals\n"
    "A.D.Becke; J.Chem.Phys.; 107, 8554, (1997)\n"
    "Implemented by Sarah Reimann ",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(b97c_en) };

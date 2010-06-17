#include "functional.h"
#include "constants.h"

template<class num>
static num new_energy(const densvars<num> &d)
{
  static const parameter DELTA = 0.1;

  num ea = d.gaa/(DELTA + pow(d.a, 4.0/3.0));
  num eb = d.gbb/(DELTA + pow(d.b, 4.0/3.0));
  return ea + eb;
}

void setup_ktx(functional &f)
{
  f.describe(XC_KTX, XC_GGA,
	     "KT exchange GGA correction",
	     "KT exchange GGA correction\n"
             "reference:\n"
             "@article{keal:3015,\n"
             "author = {Thomas W. Keal and David J. Tozer},\n"
             "collaboration = {},\n"
             "title = {The exchange-correlation potential in Kohn--Sham nuclear magnetic resonance shielding calculations},\n"
             "publisher = {AIP},\n"
             "year = {2003},\n"
             "journal = {The Journal of Chemical Physics},\n"
             "volume = {119},\n"
             "number = {6},\n"
             "pages = {3015-3024},\n"
             "keywords = {eigenvalues and eigenfunctions; ab initio calculations; bond lengths; nuclear screening; nuclear magnetic resonance; density functional theory; ionisation potential; dissociation energies},\n"
             "url = {http://link.aip.org/link/?JCP/119/3015/1},\n"
             "doi = {10.1063/1.1590634}\n"
             "}\n"
             "xcfun version: Radovan Bast (radovan.bast@uit.no)\n"
             "tested against implementation in Dalton by Dave Wilson (davidwi@kjemi.uio.no)\n"
             "compared first derivatives only\n");

  SET_GGA_ENERGY_FUNCTION(f, new_energy);

  static const double d[5] =
    {
     0.39e+02,
     0.38e+02,
     0.81e+06,
     0.82e+06,
     0.82e+06
    };

  static const double ref[21] =
    {
      0.12533312965365759737e+05, // 00000
      -.20906588852649329624e+03, // 10000
      -.22485950142470750279e+03, // 01000
      0.75553098007418890980e-02, // 00100
      0.78213561302010112253e-02, // 00010
      0.00000000000000000000e+00, // 00001
      0.12497415159317823097e+02, // 20000
      0.00000000000000000000e+00, // 11000
      -.25810603521789297204e-03, // 10100
      0.00000000000000000000e+00, // 10010
      0.00000000000000000000e+00, // 10001
      0.13794820570009042271e+02, // 02000
      0.00000000000000000000e+00, // 01100
      -.27421890417647258121e-03, // 01010
      0.00000000000000000000e+00, // 01001
      0.00000000000000000000e+00, // 00200
      0.00000000000000000000e+00, // 00110
      0.00000000000000000000e+00, // 00101
      0.00000000000000000000e+00, // 00020
      0.00000000000000000000e+00, // 00011
      0.00000000000000000000e+00  // 00002
    };

  f.add_test(XC_VARS_AB, 2, d, ref, 1e-11);
}

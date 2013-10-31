#include "xcint.hpp"

static struct alias_data aliases_array[XC_MAX_ALIASES] =
{
  {"lda","Slater exchange and VWN5 correlation",{{"slaterx",1.0},{"vwn5c",1.0}} },
  {"blyp","Becke exchange and LYP correlation",{{"beckex",1.0},{"lypc",1.0}} },
  {"pbe","PBE exchange and correlation",{{"pbex",1.0},{"pbec",1.0}} },
  {"bp86","Becke-Perdew 1986",{{"beckex",1.0},{"p86c",1.0},{"pz81c",1.0}} },
  {"kt2","Keal-Tozer 2",{{"slaterx",1.07173},{"ktx",-0.006},{"vwn5c",0.576727}} },
  {"kt3","Keal-Tozer 3",{{"slaterx",1.092},{"ktx",-0.004},{"optxcorr",-0.925452},{"lypc",0.864409}} },
  // Hybrid Functionals
  {"pbe0","Perdew-Burke-Ernzerhof",{{"pbex",0.75},{"pbec",1.0},{"exx",0.25}} }, 
  {"b3lyp","Becke-3-paramater-LYP (VWN5 form)", {{"slaterx",0.80},{"beckecorrx",0.72},{"lypc",0.81},{"vwn5c",0.19},{"exx",0.20}} }, 
  {"b3lyp-g","Becke-3-paramater-LYP (VWN3 form)", {{"slaterx",0.80},{"beckecorrx",0.72},{"lypc",0.81},{"vwn3c",0.19},{"exx",0.20}} }, 
  {"m06","M06", {{"m06c",1.0},{"m06x",1.0}} }, 
  {"m06-2x","M06 2X", {{"m06x2c",1.0},{"m06x2x",1.0}} }, 
  {"m06L","M06 L", {{"m06lc",1.0},{"m06lx",1.0}} }, 
  // Aliases for DALTON GGAKEY
  {"becke","Becke exchange",{{"beckexcorr",1.0}}},
  {"slater","Slater exchange",{{"slaterx",1.0}}},
  {"lyp","LYP correlation",{{"lypc",1.0}}},
  {"vwn","VWN5 correlation",{{"vwn5c",1.0}}},
  {"svwn","Slater exchange and VWN5 correlation",{{"slaterx",1.0},{"vwn5c",1.0}} },
  {"svwn5","Slater exchange and VWN5 correlation",{{"slaterx",1.0},{"vwn5c",1.0}} },
  {"svwn3","Slater exchange and VWN5 correlation",{{"slaterx",1.0},{"vwn3c",1.0}} },
  {"b3p86","Becke-3-paramater-LYP (VWN5 form)", {{"slaterx",0.80},{"beckecorrx",0.72},{"p86c",0.81},{"vwn5c",1.0},{"exx",0.20}} }, 
  {"b3p86-g","Becke-3-paramater-LYP (VWN5 form)", {{"slaterx",0.80},{"beckecorrx",0.72},{"p86c",0.81},{"vwn3c",1.0},{"exx",0.20}} }, 
};

struct alias_data *xcint_aliases = aliases_array;

#include "functional.hpp"

PARAMETER(XC_RANGESEP_MU) = 
{
  "Range separation inverse length [1/a0]",
  0.4
};

PARAMETER(XC_EXX) = 
{
  "Amount of exact (HF like) exchange (must be provided externally)",
  0.0
};

PARAMETER(XC_CAM_ALPHA) = 
{
  "Amount of exact (HF like) exchange within CAM-B3LYP functional",
  0.19
};

PARAMETER(XC_CAM_BETA) = 
{
  "Amount of long-range (HF like) exchange within CAM-B3LYP functional",
  0.46
};

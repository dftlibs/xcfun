#include "functional.h"
#include "constants.h"

//  Thomas-Fermi kinetic energy functional

template<class num>
static num ENERGY_FUNCTION(XC_TFK)(const densvars<num> &d)
{
    using xc_constants::CF;

    return CF*pow(d.n, 5.0/3.0);
}

NEW_LDA_FUNCTIONAL(XC_TFK);
SHORT_DESCRIPTION(XC_TFK) = "Thomas-Fermi Kinetic Energy Functional";
LONG_DESCRIPTION(XC_TFK) =
	     "Thomas-Fermi Kinetic Energy Functional\n"
	     "\n"
	     "Implemented by Andre Gomes.\n"; 
TEST_VARS(XC_TFK) = XC_A_B;
TEST_ORDER(XC_TFK) = 1;
TEST_THRESHOLD(XC_TFK) = 1e-5;
TEST_IN(XC_TFK) =  {1., .8};
TEST_OUT(XC_TFK) =   { 7.64755771625168, 7.08107195949229, 7.08107195949229, 2.62261924425641, 2.62261924425641, 2.62261924425641 };


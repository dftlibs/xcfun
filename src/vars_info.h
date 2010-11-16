#ifndef VARS_INFO_H
#define VARS_INFO_H
#include "xcfun.h"

template<int VARS> struct vars_info {};

// subset is a small set of variables contained in each set
// pot_order is the derivative order required for the potential. Set to 0 if potential is not available.

template<> struct vars_info<-1> { enum { nr_variables = 0 }; enum{ subset = -1 }; enum{ pot_order = 0 }; };
template<> struct vars_info<XC_A> { enum { nr_variables = 1 }; enum{ subset = -1 }; enum{ pot_order = 1 };};
template<> struct vars_info<XC_N> { enum { nr_variables = 1 }; enum{ subset = -1 }; enum{ pot_order = 1 };};
template<> struct vars_info<XC_A_B>   { enum { nr_variables = 2 }; enum{ subset = XC_A }; enum{ pot_order = 1 };};
template<> struct vars_info<XC_N_S>   { enum { nr_variables = 2 }; enum{ subset = XC_N }; enum{ pot_order = 1 };};
template<> struct vars_info<XC_A_GAA> { enum { nr_variables = 2 }; enum{ subset = XC_A }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_N_GNN> { enum { nr_variables = 2 }; enum{ subset = XC_N }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_A_GAA_LAPA>  { enum { nr_variables = 3 }; enum{ subset = XC_A_GAA }; enum{ pot_order = 2 };};
template<> struct vars_info<XC_A_GAA_TAUA> { enum { nr_variables = 3 }; enum{ subset = XC_A_GAA }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_N_GNN_LAPN>  { enum { nr_variables = 3 }; enum{ subset = XC_N_GNN }; enum{ pot_order = 2 };};
template<> struct vars_info<XC_N_GNN_TAUN> { enum { nr_variables = 3 }; enum{ subset = XC_N_GNN }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_A_B_GAA_GAB_GBB> { enum { nr_variables = 5 }; enum{ subset = XC_A_B }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_N_S_GNN_GNS_GSS> { enum { nr_variables = 5 }; enum{ subset = XC_N_S }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_A_B_GAA_GAB_GBB_LAPA_LAPB> { enum { nr_variables = 7 }; enum{ subset = XC_A_B_GAA_GAB_GBB }; enum{ pot_order = 2 };};
template<> struct vars_info<XC_A_B_GAA_GAB_GBB_TAUA_TAUB> { enum { nr_variables = 7 }; enum{ subset = XC_A_B_GAA_GAB_GBB }; enum{ pot_order = 0 };};
template<> struct vars_info<XC_N_S_GNN_GNS_GSS_LAPN_LAPS> { enum { nr_variables = 7 }; enum{ subset = XC_N_S_GNN_GNS_GSS }; enum{ pot_order = 2 };};
template<> struct vars_info<XC_N_S_GNN_GNS_GSS_TAUN_TAUS> { enum { nr_variables = 7 }; enum{ subset = XC_N_S_GNN_GNS_GSS }; enum{ pot_order = 0};};
/*template<> struct vars_info<XC_A_B_GAX_GAY_GAZ_GBX_GBY_GBZ> { enum { nr_variables = 9 }; enum{ subset = XC_A_B }; };
template<> struct vars_info<XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB> { enum { nr_variables = 9 }; enum{ subset = XC_A_B }; };
template<> struct vars_info<XC_A_B_GAX_GAY_GAZ_GBX_GBY_GBZ_LAPA_LAPB> { enum { nr_variables = 10 }; enum{ subset = XC_A_B }; };
template<> struct vars_info<XC_A_B_GAX_GAY_GAZ_GBX_GBY_GBZ_TAUA_TAUB> { enum { nr_variables = 10 }; enum{ subset = XC_A_B }; };
template<> struct vars_info<XC_A_B_GAX_GAY_GAZ_GBX_GBY_GBZ_LAPA_LAPB_TAUA_TAUB> { enum { nr_variables = 12 }; 
enum{ subset = XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB }; }; */


#endif

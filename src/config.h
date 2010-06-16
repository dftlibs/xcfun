#ifndef CONFIG_H
#define CONFIG_H

//Maximum derivative order. Lower orders can be generated
//for GGA's and MGGA's, to avoid huge code size.
#define XC_MAX_ORDER 2
#define XC_LDA_MAX_ORDER XC_MAX_ORDER
#define XC_GGA_MAX_ORDER XC_MAX_ORDER
#define XC_MGGA_MAX_ORDER XC_MAX_ORDER

#endif

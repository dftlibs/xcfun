#ifndef XCFUN_H
#define XCFUN_H

#define XCFUN_API_VERSION 0

// Which set (mode) of variables to use
#define XC_VARS_A   0 // 100% spin polarized, i.e. only alpha variables
#define XC_VARS_N   1 // No spin polarization, total density variables
#define XC_VARS_AB  2 // Alpha/beta variables
#define XC_VARS_NS  3 // Total density/spin density variables
#define XC_NR_MODES 4

// Type of functional
#define XC_LDA      0 // Local density
#define XC_GGA      1 // Local density & gradient
#define XC_MGGA     2 // Local density, gradient and kinetic energy density
#define XC_MLGGA    3 // Local density, gradient, laplacian and kinetic energy density
#define XC_NR_TYPES 4

#ifdef __cplusplus
extern "C" {
#endif

  double xcfun_version(void);
  const char *xcfun_splash(void);
  int xcfun_test(void);
  
  typedef struct xc_functional_data * xc_functional;

  xc_functional xc_new_functional(void);
  void xc_free_functional(xc_functional fun);

  // Evaluate the functional at density
  void xc_eval(xc_functional fun, int order,
	       const double *density,
	       double *result);
  /* Vector version of xc_eval. 
     density_pitch = density[start_of_second_point] - density[start_of_first_point],
     likewise for result_pitch.
   */
  void xc_eval_vec(xc_functional fun, int order, int nr_points,
		   const double *density,
		   int density_pitch,
		   double *result,
		   int result_pitch);

  // Calculate the xc potential for fun. This is currently only supported for some
  // types of functionals (LDA or GGA, AB mode). Here density should be the normal
  // density as for xc_eval, appended with the extra data needed to construct the
  // potential. For GGA's this is the laplacian of the density. The energy density
  // is computed as a by-product.
  void xc_potential(xc_functional fun, const double *density, double *e_xc, double *v_xc);

  // Which variables to use/differentiatiate wrt to
  void xc_set_mode(xc_functional fun, int mode);

  // The type of the currently defined functional (LDA, GGA etc)
  int xc_get_type(xc_functional fun);

  // The highest order supported by the currently defined functional.
  // This depends on compile time parameters, it's in principle unlimited.
  int xc_max_order(xc_functional fun);

  // Length of the density[] argument to eval()
  int xc_input_length(xc_functional fun);

  // Length of the result[] argument to eval()
  int xc_output_length(xc_functional fun, int order);

  // Index into result[] for derivative with given index (length as input_length() )
  int xc_derivative_index(xc_functional fun, const int derivative[]);

  // List of all settings, there are XC_NR_PARAMS settings in total.
#ifndef XCFUN_INTERNAL
#include "xcfun_autogen.h"
#endif

  /* Discover and manipulate settings */

  // Return the internal name of setting param
  const char *xc_name(int param);
  // Describe in one line what the setting does
  const char *xc_short_description(int param);
  // Long description of the setting, ends with a \n
  const char *xc_long_description(int param);
  // Is this setting a functional?
  int xc_is_functional(int param);
  // Set the setting
  void xc_set_param(xc_functional fun, int param, double value);
  // Get the current value of the setting.
  double xc_get_param(xc_functional fun, int param);

  /* Transform the output of xc_eval to a different mode,
     for example because your program wants AB mode but you
     want to take advantage of knowing that you only care
     about non-polarizing N mode derivatives. */
  void xc_transform(int order,
		    int from_mode, const double *from_data,
		    int to_mode, double *to_data);

#ifdef __cplusplus
} // End of extern "C"
#endif

// Derivative indices into xc_eval output

#define XC_D0 0
#define XC_D1 1
#define XC_D2 2

#define XC_D00 0
#define XC_D10 1
#define XC_D01 2
#define XC_D20 3
#define XC_D11 4
#define XC_D02 5

#define XC_D00000 0
#define XC_D10000 1
#define XC_D01000 2
#define XC_D00100 3
#define XC_D00010 4
#define XC_D00001 5
#define XC_D20000 6
#define XC_D11000 7
#define XC_D10100 8
#define XC_D10010 9
#define XC_D10001 10
#define XC_D02000 11
#define XC_D01100 12
#define XC_D01010 13
#define XC_D01001 14
#define XC_D00200 15
#define XC_D00110 16
#define XC_D00101 17
#define XC_D00020 18
#define XC_D00011 19
#define XC_D00002 20

#endif

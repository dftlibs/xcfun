#ifndef XCFUN_H
#define XCFUN_H

#define XCFUN_API_VERSION 0

// Which set (mode) of variables to use
#define XC_VARS_A   0 // 100% spin polarized, i.e. only alpha variables
#define XC_VARS_R   1 // No spin polarization, total density variables
#define XC_VARS_AB  2 // Alpha/beta variables
#define XC_VARS_RS  3 // Total density/spin density variables
#define XC_NR_MODES 4

// Type of functional
#define XC_LDA      0 // Local density
#define XC_GGA      1 // Local density & gradient
#define XC_MGGA     2 // Local density, gradient and kinetic energy density
//#define XC_M2GGA    3 // Local density, gradient, laplacian and kinetic energy density
#define XC_NR_TYPES 3

double xcfun_version(void);
const char *xcfun_splash(void);
int xcfun_test(void);

class xc_functional
{
 public:
  xc_functional(void);
  ~xc_functional(void);

  // If density would result in infinite derivatives,
  // make a tiny modification to density in order to 
  // get finite values. The modified density is expected
  // to give rise to the same derivative values for those
  // derivatives that are not too large, and finite but very
  // large derivatives in the case where the unregularized
  // density would give infinities.
  void regularize_density(double *density);

  // Evaluate the functional at density
  void eval(double *result, int order, const double *density);

  // Which variables to use/differentiatiate wrt to
  void set_mode(int mode);

  // The type of the currently defined functional (LDA, GGA etc)
  int get_type(void) const;

  // The highest order supported by the currently defined functional.
  // This depends on compile time parameters, it's in principle unlimited.
  int get_max_order(void) const;

  // Length of the density[] argument to eval()
  int input_length(void) const;

  // Length of the result[] argument to eval()
  int output_length(int order) const;

  // Index into result[] for derivative with given index (length as input_length() )
  int derivative_index(const int derivative[]) const;

  // Return the name of setting n, or NULL if there is no such setting.
  // Settings are numbered consecutively from 0
  const char *setting_name(int n) const;
  // Set name = value, return 0 if name is a valid setting
  int set_setting(const char *name, double value);
  double get_setting(const char *name) const;
  bool is_set(const char *name) const;
  bool is_functional(const char *name) const;
  const char *setting_short_description(const char *name) const;
  const char *setting_long_description(const char *name) const;
  class xc_functional_data;
 protected:
  xc_functional_data *d;
};

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

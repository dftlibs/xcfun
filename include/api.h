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
#define XC_NR_TYPES 3

const char *xcfun_reference(void);

class xc_functional
{
 public:
  void set_mode(int mode);
  int set_parameter(const char *parameter, double value);
  void eval(double *result, int order, const double *density);
  void regularize_density(double *density);

  int get_type(void) const;
  void set_type(int type);
  int get_max_order(void) const;

  int input_length(void) const;
  int output_length(int order) const;
  int derivative_index(const int derivative[]) const;

  int nr_parameters(void) const;
  bool is_a_functional(const char *parameter);
  const char *parameter_name(int parameter_nr) const;
  const char *describe_parameter(const char *parameter) const;
  const char *describe_implementation(const char *parameter) const;
};

#endif

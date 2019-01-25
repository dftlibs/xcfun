%module xcfun_swig

%{
#define SWIG_FILE_WITH_INIT
#include "XCFun/xcfun.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%inline %{
int xc_eval_setup_swig(xc_functional fun, enum xc_vars vars, enum xc_mode mode, int order)
    {
       return xc_eval_setup(fun, vars, mode, order);
    }
%}

%apply (int DIM1, double* IN_ARRAY1 ) {(int size_density, const double* density)}
%apply (int DIM1, double* INPLACE_ARRAY1 ) {(int size_result, double* result)}

%inline %{
void xc_eval_swig(xc_functional fun, int size_density, const double *density, int size_result, double *result)
    {
        xc_eval(fun, density, result);
    }
%}

%apply (int DIM1, int DIM2, double* IN_ARRAY2 ) {(int nr_points, int size_density, const double *density)}
%apply (int DIM1, int DIM2, double* INPLACE_ARRAY2 ) {(int nr_points_res, int size_result, double *result)}

%inline %{
void xc_eval_vec_swig(xc_functional fun, int nr_points, int size_density, const double *density,
                                         int nr_points_res, int size_result, double *result)
    {
        xc_eval_vec(fun, nr_points, density, size_density, result, size_result);
    }
%}

%inline %{
double xcfun_version_swig(void)
    {
       return xcfun_version();
    }
%}

%inline %{
const char* xcfun_splash_swig(void)
    {
       return xcfun_splash();
    }
%}

%ignore xc_deserialize;
%ignore xc_derivative_index;
%ignore xc_set_fromstring;

%ignore xc_eval_setup;
%ignore xc_eval;
%ignore xc_eval_vec;

%ignore xc_get;
%ignore xc_describe_long;
%ignore xc_serialize;
%ignore xc_deserialize;

%include "xcfun.h"

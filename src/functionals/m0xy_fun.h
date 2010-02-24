#ifndef M0xy_H
#define M0xy_H

// common functions for MO5 and M06 family of (hybrid) meta-gga functionals
//

namespace m0xy_metagga_xc_internal
{

  const char *m05_2x_x_reference = 
    "M05-2X Meta-GGA Exchange Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m05_x_reference = 
    "M05 Meta-GGA Exchange Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m06_l_x_reference = 
    "M06-L Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "J. Chem. Phys. 125, 194101 (2006)\n"; 

  const char *m06_2x_x_reference = 
    "M06-2X Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 

  const char *m06_x_reference = 
    "M06 Meta-GGA Exchange Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 


  const char *m05_2x_c_reference = 
    "M05-2X Meta-GGA Correlation Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m05_c_reference = 
    "M05 Meta-GGA Correlation Functional\n"
      "Y Zhao, N. E. Schultz and D. G. Truhlar, "
      "J. Chem. Theory Comput. 2, 364 (2006)\n"; 
 
  const char *m06_l_c_reference = 
    "M06-L Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "J. Chem. Phys. 125, 194101 (2006)\n"; 

  const char *m06_2x_c_reference = 
    "M06-2X Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 

  const char *m06_c_reference = 
    "M06 Meta-GGA Correlation Functional\n"
      "Y Zhao and D. G. Truhlar, "
      "Theor. Chem. Account 120, 215 (2008)\n"; 


// additional information from the authors, and reference implementations 
// in F77, can be found at
//
//    http://comp.chem.umn.edu/mfm/
//


// the alpha below is defined at: 

   static num alpha_m0xc = 1.0;

// rho is the density, tau the kinetic energy density, see TCA 120 215 (2008), page 219 
// zet : eq (3) of paper above
  template<class num>
  static num zet(const num &rho, const num &tau)
  {
    using xc_constants::CF;

    return ((2.0*tau/pow(rho,5.0/3.0)) - CF);
  }

  template<class num>
  static num gamma(const num &chi, const num &zet)
  {
    return 1 + alpha_m0xc*(chi*chi + zet); 
  }

  template<class num>
  static num h(const static double d, const num &chi, const num &zet)
  {
    num gam1 = gamma(chi,zet);

    num t1 =  d[0]/gam1;
    num t2 = (d[1]*pow(chi,2.0) + d[2]*zet)/pow(gam1, 2.0);
    num t3 = (d[3]*pow(chi,4.0) + d[4]*pow(chi,2.0)+ d[5]*pow(zet, 2.0))/pow(gam1, 3.0);

    return t1+t2+t3;
  }

  template<class num>
  static num fw(const static double a, const num &rho, const num &tau)
  {
    using pw91_like_x_internal::pw91k_prefactor;

    num tau_lsda = pw91k_prefactor(rho);

    num  t = tau_lsda/tau;
    num  w = (t - 1)/(t + 1); 
    num fw = 0.0;

    for (int i = 0; i <= 11; i++)
       fw += a[i] * pow(w,i);

    return fw;
  }

}

#endif

// This file is ment to be included from xc_fun.cpp 

/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<class num>
static num lyp_correlation(const num &R,
			   const num &rhoa,
			   const num &rhob,
			   const num &Z,
			   const num &grada2,
			   const num &gradb2)
{
    const double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    using xc_constants::CF;

    num rho2 = R*R;
    num rhom13 = pow(R,-1.0/3.0);
    num denom = 1+D*rhom13;
    num omega = exp(-C*rhom13)/denom*pow(R,-11.0/3.0);
    num delta = rhom13*(C + D/denom);

    num t1 = pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +
				   pow(rhob,8.0/3.0));
    t1 +=  ((47.0 - 7.0*delta)/18.0)*Z;
    t1 += -(2.5 -delta/18.0)*(grada2+gradb2);
    t1 +=  (11.0-delta)/9.0*(rhoa*grada2 + rhob*gradb2)/R;
    num t5 = -2.0/3.0*rho2*Z;
    t5 += ((2.0/3.0*rho2-rhoa*rhoa)*gradb2 +
	   (2.0/3.0*rho2-rhob*rhob)*grada2);
    return -A*(4*rhoa*rhob/(denom*R)
	      +B*omega*(rhoa*rhob*t1+t5));
}


template<class num>
static num becke_exchange(const num &rhoa,
			  const num &rhob,
			  const num &grada2,
			  const num &gradb2)
{
  static const double BECKE_THRESHOLD = 1e-14;
  static const double BETA = 0.0042;  
  num ea = 0, eb = 0;

  if (rhob>BECKE_THRESHOLD)
    {
      num xb = sqrt(gradb2)*pow(rhob,-4.0/3.0);
      num rb = pow(rhob,4.0/3.0);
      num denomb = 1.0 +6.0*xb*BETA*asinh(xb);
      eb = rb*xb*xb/denomb; 
    } 
  if (rhoa>BECKE_THRESHOLD) 
    {
      num xa = sqrt(grada2)*pow(rhoa,-4.0/3.0);
      num ra = pow(rhoa,4.0/3.0);
      num denoma = 1.0 +6.0*BETA*xa*asinh(xa);
      ea = ra*xa*xa/denoma;
    }
  return -BETA*(ea+eb);
}

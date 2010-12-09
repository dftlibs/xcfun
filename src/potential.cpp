#include "xcfun_internal.h"

// Compute the xc potential(s). This is simple for lda and more complicated for 
// GGA's. For metaGAA's that depend on tau this is a non-local problem and 
// cannot be solved by xcfun. Laplacian dependent metaGGA's could be implemented..
// Another option is to evaluate the divergence directly in cartesian coordinates,
// but one have to find a way to do this without introducing the gradient components
// explicitly (for niceness)
void xc_potential(xc_functional fun, const double *density, double *e_xc, double *v_xc)
{
  if (fun->mode == XC_VARS_AB)
    {
      if (fun->type == XC_LDA) 
	{
	  double out[3];
	  xc_eval(fun,1,density,out);
	  e_xc[0] = out[0];
	  v_xc[0] = out[1];
	  v_xc[1] = out[2];
	}
      // Expecting lap_a and lap_b as element 5 and 6 of density
      // Using that v_xc = dE/dn - div dE/dgradn
      // When deriving the expressions it helps to have the operator
      // (g_a).nabla act on the basic variables (and the same for g_b). 
      else if (fun->type == XC_GGA)  
	{
	  const int gaa = 2, gab = 3, gbb = 4, lapa = 5, lapb = 6;
	  double out[21];
	  xc_eval(fun,2,density,out);
	  e_xc[0] = out[XC_D00000];

	  v_xc[0]  = out[XC_D10000];
	  v_xc[0] -= 2*density[lapa]*out[XC_D00100] + density[lapb]*out[XC_D00010];
	  v_xc[0] -= 2*(out[XC_D10100]*density[gaa]   + 
			out[XC_D01100]*density[gab] +
			out[XC_D00200]*(2*density[lapa]*density[gaa]) +
			out[XC_D00110]*(density[lapa]*density[gab] + density[lapb]*density[gaa]) +
			out[XC_D00101]*(2*density[lapb]*density[gab]) 
			); 
	  v_xc[0] -= (out[XC_D10010]*density[gab] +
		      out[XC_D01010]*density[gbb] +
		      out[XC_D00110]*(2*density[lapa]*density[gab]) +
		      out[XC_D00020]*(density[lapb]*density[gab] + density[lapa]*density[gbb]) +
		      out[XC_D00011]*(2*density[lapb]*density[gbb])); 

	  v_xc[1]  = out[XC_D01000];
	  v_xc[1] -= 2*density[lapb]*out[XC_D00001] + density[lapa]*out[XC_D00010];
	  v_xc[1] -= 2*(out[XC_D01001]*density[gbb]   + 
			out[XC_D10001]*density[gab]  +
			out[XC_D00002]*(2*density[lapb]*density[gbb]) +
			out[XC_D00011]*(density[lapb]*density[gab] + density[lapa]*density[gbb]) +
			out[XC_D00101]*(2*density[lapa]*density[gab])  ); 
	  v_xc[1] -= (out[XC_D01010]*density[gab] +
		      out[XC_D10010]*density[gaa] +
		      out[XC_D00011]*(2*density[lapb]*density[gab]) +
		      out[XC_D00020]*(density[lapa]*density[gab] + density[lapb]*density[gaa]) +
		      out[XC_D00110]*(2*density[lapa]*density[gaa])); 
	}
      else
	{
	  xcint_die("xc_potential() not implemented for metaGGA's",fun->type);
	}
    }
  else
    {
      if (fun->type == XC_LDA) 
	{
	  double out[2];
	  xc_eval(fun,1,density,out);
	  e_xc[0] = out[0];
	  v_xc[0] = out[1];
	}      
      else
	{
	  xcint_die("xc_potential() GGA only implemented for AB mode",
		    fun->mode);
	}
    }
}


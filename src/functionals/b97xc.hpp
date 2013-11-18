#ifndef B97XC_H
#define B97XC_H

namespace b97xc
{
    
    
    // square of spin-density gradient
    template<class num>
    static num spin_dens_gradient_ab2(const num &gaa, const num &a_43){
        
        return abs(gaa)/a_43/a_43;
        
    }
    
   
    template<class num>
    static num ux_ab(const parameter &Gamma, const num &spin_dens_grad){
        
        
        return Gamma*spin_dens_grad/(1+Gamma*spin_dens_grad);
    }
    
    template<class num>
    static num enhancement(const parameter &Gamma, const parameter c_params[], const num &spin_dens_grad)
    {
        
        num ux = ux_ab(Gamma, spin_dens_grad);
        
        return  c_params[0] + c_params[1]*ux + c_params[2]*ux*ux;
        
    }
    
    
}


#endif

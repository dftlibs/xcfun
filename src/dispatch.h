


// _adds_ the functional results weighted with weight.
typedef void (*wrapped_functional)(void *result, const void *dv, const void *weight);

// Dispatch table for all supported functional types and orders
static wrapped_functional
double_tab[XC_FUNLIST_LEN][XC_MODES_LEN][XC_MAX_ORDER+1] = {{{0}}};
static const char *fun_name_tab[XC_FUNLIST_LEN] = {0}; 


// This is a triple loop, over all functionals, modes
// and degrees. It tries carefully not to mention
// functional::energy with parameters that are not
// needed, in order to generate as little code as possible.

template<int fun_id,  int mode_mask, int Ndeg>
struct table
{
  template<class scalar>
  static void wrapped_call(void *result, const void *dv, const void *weight)
  {
    typedef taylor<scalar,xc_mode_info<lowest_set_bit<mode_mask>::bit>::Nvar,Ndeg> ttype;
    const densvars<ttype> *d = reinterpret_cast<const densvars<ttype> *>(dv); 
    ttype *res = reinterpret_cast<ttype *>(result);
    *res += *reinterpret_cast<const scalar *>(weight)*functional<fun_id>::energy(*d);
  }
  static void setup(void)
  {
    double_tab[fun_id][bitpos<lowest_set_bit<mode_mask>::bit>::pos][Ndeg] = 
      table<fun_id,lowest_set_bit<mode_mask>::bit,Ndeg>
      ::template wrapped_call<double>;
    table<fun_id,mode_mask,Ndeg-1>::setup();
  }
};

// End of order loop
template<int fun_id, int mode_mask>
struct table<fun_id,mode_mask,-1>
{
  static void setup(void) 
  {
    table<fun_id,mode_mask - lowest_set_bit<mode_mask>::bit, 
      functional<fun_id>::max_order>::setup();
  }
};

// End of mode loop
template<int fun_id, int Ndeg>
struct table<fun_id,0,Ndeg>
{
  static void setup(void) 
  { 
    fun_name_tab[fun_id] = functional<fun_id>::get_name();
    table<fun_id+1, functional<fun_id+1>::supported_modes,
      functional<fun_id+1>::max_order>::setup();
  }
};

// Dummy functional to terminate loops
template<>
struct functional<XC_FUNLIST_LEN>
{
  enum { supported_modes = 0 };
  enum { max_order = 0 };
};

// End of functional loop
template<>
struct table<XC_FUNLIST_LEN,0,0>
{
  static void setup(void) {}
};

void make_sure_tables_are_set_up(void)
{
  static int is_setup = 0;
  if (!is_setup)
    {
      table<0,functional<0>::supported_modes,
	functional<0>::max_order>::setup();
      is_setup = 1;
    }
}

struct functional_term
{
  int fun_id;
  double weight;
  functional_term *pnext;
}


void xc_set_functional(enum xc_mode mode, 
		       const double params[XC_PARAMS_LEN])
{
  make_sure_tables_are_set_up();
  current_mode = mode;
  for (int i=0;i<XC_PARAMS_LEN;i++)
    xc_params[i] = params[i];
}

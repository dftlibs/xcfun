#ifndef TEMPLATE_MISC_H
#define TEMPLATE_MISC_H

template<int N>
struct lowest_set_bit
{
  enum { bit = (N & 1) ? 1 : 2*lowest_set_bit<N/2>::bit };
};

template<>
struct lowest_set_bit<0>
{
  enum { bit = 0 };
};

template<int N>
struct bitpos
{
  enum { pos = (N & 1) ? 0 : 1+bitpos<N/2>::pos };
};

template<>
struct bitpos<0>
{
  enum { pos = -1000 };
};
    
#endif

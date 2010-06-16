/* A c++ program that can be linked without libstdc++.
Requires compilation with -fno-rtti -fno-exceptions
with gcc and intel, similar options with other compilers?
 */

#include <cstdio> //Seems to be ok, requires libc
#include <cstdlib>
#include "struct_array.h"

class C
{
public:
  C(int x) { data = x; mem = (int *)malloc(sizeof*mem); }
  ~C() { printf("Destroying C object with data = %i\n",data); free(mem); }
protected:
  int data;
  int *mem;
};

//Templates are ok
template<int N>
int retN()
{
  // This is ok (normal static data)
  static double f[] = {1,2,3};
  // Creating a static object bring in locking
  /* Brings in __cxa_guard_acquire/release
     static C cintemplate(19);  */
  return N+f[0];
}

static C cstatic(5); //This is also ok with g++/gcc

int main(void)
{
  C cstack(12);
  array<int> a;
  a.construct();
  a.push_back(12);
  printf("from array: %i\n",a[0]);
  a[0] = 18;
  printf("from array: %i\n",a[0]);
  a.destruct();
  /* Requires libstdc++
  C *cp = new C(17);
  delete cp;
  */
  printf("ret12: %i\n",retN<12>());
  return 0;
}

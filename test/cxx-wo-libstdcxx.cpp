/* A c++ program that can be linked without libstdc++ */

#include <cstdio> //Seems to be ok, requires libc
#include <cstdlib>

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
  /* Requires libstdc++
  C *cp = new C(17);
  delete cp;
  */
  printf("ret12: %i\n",retN<12>());
  return 0;
}

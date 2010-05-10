#ifndef ARRAY_H
#define ARRAY_H

// Stupid array class because we don't want to rely on libstdc++

#include <cassert>

template<typename T>
class array
{
public:
  array(void)
  {
    nr_items = 0;
    max_items = 2;
    items = new T[max_items];
  }
  array(const array &a)
  {
    nr_items = a.nr_items;
    max_items = a.max_items;
    items = new T[max_items];
    for (int i=0;i<nr_items;i++)
      items[i] = a[i];
  }
  array<T> &operator=(const array &a)
  {
    if (this != &a)
      {
	if (a.nr_items>max_items)
	  grow(a.nr_items);
	nr_items = a.nr_items;
	for (int i=0;i<nr_items;i++)
	  items[i] = a[i];
      }
    return *this;
  }
  ~array(void)
  {
    delete[] items;
  }
  int size(void) const
  {
    return nr_items;
  }
  T &operator[](int i)
  {
    assert(i>=0);
    assert(i<nr_items);
    return items[i];
  }
  const T &operator[](int i) const
  {
    assert(i>=0);
    assert(i<nr_items);
    return items[i];
  }
  void push_back(const T &x)
  {
    grow(2*max_items);
    items[nr_items++] = x;
  }
  int find(const T &x)
  {
    for (int i=0;i<nr_items;i++)
      if (items[i] == x)
	return i;
    return -1;
  }
  void clear(void)
  {
    delete[] items;
    nr_items = 0;
    max_items = 2;
    items = new T[max_items];
  }
  void extend(int len, const T &fill)
  {
    grow(len);
    for (int i=nr_items;i<len;i++)
      items[i] = fill;
    nr_items = len;
  }
protected:
  // Allocated for at least capacity items
  void grow(int capacity)
  {
    if (capacity > max_items)
      {
	max_items = capacity;
	T *new_items = new T[max_items];
	for (int i=0;i<nr_items;i++)
	  new_items[i] = items[i];
	delete[] items;
	items = new_items;
      }
  }
  int nr_items;
  int max_items;
  T *items;
};


#endif

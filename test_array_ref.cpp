/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_array_ref test_array_ref.cpp
*/

#include <array_ref>

void foo(int A[] , size_t N); // Traditional API
void foo(const int A[], size_t N); // Traditional API

void foo(array_ref<int[]> A); // Reference API
void foo(array_ref<const int[]> A); // Reference API

void
bar()
{
  enum { L = ... };
  int buffer[L];
  array_ref<int[]> A(buffer, L);

  assert(L == A.size());
  assert(& A[0] == buffer);

  foo(array);
}

int
main()
{
}


/*

<http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43622>
This test-case for GCC bug 43622 shows how
typeinfo for non-complex 128-bit types DOES NOT work, but
typeinfo for complex 128-bit types DOES work.

$ g++-4.5 --version
g++-4.5 (Ubuntu/Linaro 4.5.1-7ubuntu2) 4.5.1
Copyright (C) 2010 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

$ uname --all
Linux dev8 2.6.35-28-generic #50-Ubuntu SMP Fri Mar 18 18:42:20 UTC 2011 x86_64 GNU/Linux

$ g++-4.5 gcc_bug_43622.cpp && ./a.out
complex_int128 = Cn
complex_uint128 = Co
complex_binary128 = Cg

*/

typedef int int128 __attribute__((mode(TI)));
typedef unsigned int uint128 __attribute__((mode(TI)));
typedef float binary128 __attribute__((mode(TF)));

typedef _Complex int complex_int128 __attribute__((mode(CTI)));
typedef _Complex unsigned int complex_uint128 __attribute__((mode(CTI)));
typedef _Complex float complex_binary128 __attribute__((mode(TC)));

#include <iostream>

#include <typeinfo>

using namespace std;

int main()
{
/* Fixed:
/tmp/ccBfQl0L.o: In function `main':
gcc_bug_43622.cpp:(.text+0xa): undefined reference to `typeinfo for __int128'
gcc_bug_43622.cpp:(.text+0x41): undefined reference to `typeinfo for unsigned __int128'
gcc_bug_43622.cpp:(.text+0x78): undefined reference to `typeinfo for __float128'
collect2: ld returned 1 exit status
*/
        cout << '\n';
	cout << "int128 = " << typeid(int128).name() << '\n';
	cout << "uint128 = " << typeid(uint128).name() << '\n';
	cout << "binary128 = " << typeid(binary128).name() << '\n';

        cout << '\n';
	cout << "complex_int128 = " << typeid(complex_int128).name() << '\n';
	cout << "complex_uint128 = " << typeid(complex_uint128).name() << '\n';
	cout << "complex_binary128 = " << typeid(complex_binary128).name() << '\n';

        cout << '\n';
	cout << "__int128 = " << typeid(__int128).name() << '\n';
	cout << "unsigned __int128 = " << typeid(unsigned __int128).name() << '\n';
	cout << "__float128 = " << typeid(__float128).name() << '\n';

        cout << '\n';
	cout << "_Complex __int128 = " << typeid(_Complex __int128).name() << '\n';
	cout << "_Complex unsigned __int128 = " << typeid(_Complex unsigned __int128).name() << '\n';
	cout << "_Complex __float128 = " << typeid(_Complex __float128).name() << '\n';

	return 0;
}

#include <spheroidal_harmonic.h>

#include <iostream>
#include <iomanip>

void
test_oblate_harmonic()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();

  int mode = 0;
  int nuevo = 5;
  double R[6], T[6];
  for (int l = 0; l <= 5; ++l)
    {
      std::cout << " l = " << l << '\n';
      for (int m = 0; m <= 5; ++m)
	{
	  std::cout << " m = " << m << '\n';
	  for (int i = 1; i <= 20; ++i)
	    {
	      double x = i * 0.1;
	      obl_sph_harm(x, m, l, mode, R, T, nuevo);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << R[l]
			<< ' ' << std::setw(w) << T[l] << '\n';
	    }
	}
    }
}

void
test_prolate_harmonic()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();

  int mode = 0;
  int nuevo = 5;
  double P[6], Q[6];
  for (int l = 0; l <= 5; ++l)
    {
      std::cout << " l = " << l << '\n';
      for (int m = 0; m <= 5; ++m)
	{
	  std::cout << " m = " << m << '\n';
	  for (int i = 1; i <= 20; ++i)
	    {
	      double x = 1.0 + i * 0.1;
	      pro_sph_harm(x, m, l, mode, P, Q, nuevo);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << P[l]
			<< ' ' << std::setw(w) << Q[l] << '\n';
	    }
	}
    }
}

int
main()
{
  std::cout << "\n Oblate Harmonics\n";
  test_oblate_harmonic();

  std::cout << "\n Prolate Harmonics\n";
  test_prolate_harmonic();
}

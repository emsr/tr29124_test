/**
 * A function to print all prime factors of a given number n
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/sf_prime.h>

void
print_prime_factors(unsigned int n)
{
  // Print the number of 2s that divide n
  while (n % 2 == 0)
    {
      std::cout << ' ' << 2;
      n /= 2;
    }
  const int n_max = std::sqrt(n);
  // n must be odd at this point.  So we can skip
  // one element (Note i = i +2)
  for (int i = 3; i <= n_max; i += 2)
    {
      // While i divides n, print i and divide n
      while (n % i == 0)
        {
          std::cout << ' ' << i;
          n /= i;
        }
    }

  // This condition is to handle the case when n
  // is a prime number greater than 2
  if (n > 2)
    std::cout << ' ' << n;
  std::cout << '\n';
}

void
prime_factors(unsigned int n)
{
  const unsigned long p_max = std::sqrt(n);
  unsigned long i = 0;
  auto p = emsr::prime(i);
  while (p > 0 && p <= p_max)
    {
      while (n % p == 0)
        {
          std::cout << ' ' << p;
          n /= p;
        }
      p = emsr::prime(++i);
    }
  if (n > 1)
    std::cout << ' ' << n;
  std::cout << '\n';
}

// Write 16-bit compressed primes.
void
write_primes()
{
  const auto max = std::numeric_limits<unsigned short>::max();
  unsigned short index = 0;
  for (unsigned int i = 0; i < emsr::num_primes(); ++i)
    {
      auto p = emsr::prime(i);
      if (p >= max)
	{
	  if (index == 0)
	    index = i;
	  p -= max;
	}
      std::cout << std::setw(7) << p
		<< "u,";
      if ((i + 1) % 8 == 0)
        std::cout << '\n';
    }
  std::cout << index << '\n';
}

int
main()
{
  while (true)
    {
      std::cout << "\nEnter integer: ";
      unsigned int n;
      std::cout << '\n';
      std::cin >> n;
      if (std::cin.bad() || std::cin.fail())
	break;
      print_prime_factors(n);
      prime_factors(n);
    }

  //write_primes();
}

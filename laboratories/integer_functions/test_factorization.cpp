/**
 *
 */

#include <iostream>
#include <cmath>

// Fermat's method for factoring. Based on Algorithm C in
// @see Donald E. Knuth "The Art of Computer Programming", vol. 2, 3rd ed,
// Seminumerical Algorithms, Addison-Wesley, pp.386-388
std::pair<unsigned int, unsigned int>
factor_fermat(unsigned int N)
{
  unsigned int sq = int(std::sqrt(N));
  unsigned int xp = 2 * sq + 1;
  unsigned int yp = 1;
  int r = sq * sq - N;

  do
    {
      while (r > 0)
	{
	  r -= yp;
	  yp += 2;
	} 
      while (r < 0)
	{
	  r += xp;
	  xp += 2;
	}
    }
  while (r != 0);

  return std::make_pair((xp - yp) / 2, (xp + yp - 2) / 2);
}

int
main()
{
  unsigned int N;
  std::cout << "\nInput a positive, odd integer: ";
  std::cin >> N;
  if (std::cin.bad() || std::cin.fail())
    return 1;

  auto facts = factor_fermat(N);

  std::cout << '\n' << ' ' << facts.first << ' ' << facts.second << '\n';

  return 0;
}

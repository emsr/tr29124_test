
#include <iostream>

#include <test/imag/test_imag.h>
#include <test/real/test_real.h>
#include <test/spot/test_spot.h>
#include <utility/util_timer.h>

namespace
{
  static void test_real_imag(void)
  {
    const Util::timer tm;
    const bool test_real_ok = test::real::test_real(true);
    const bool test_imag_ok = test::imag::test_imag(true);
    const double elapsed = tm.elapsed();

    if(test_real_ok)
    {
      std::cout << "Real test: Pass: All tests OK." << std::endl;
    }
    else
    {
      std::cout << "Real test: Fail: Not all tests OK." << std::endl;
    }

    if(test_imag_ok)
    {
      std::cout << "Imag test: Pass: All tests OK." << std::endl;
    }
    else
    {
      std::cout << "Imag test: Fail: Not all tests OK." << std::endl;
    }

    std::cout << "Elapsed time: " << elapsed << " (seconds)" << std::endl;
  }
}

int main(int, char**)
{
//  test::spot::test_spot();

  ::test_real_imag();
}

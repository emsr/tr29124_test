#include <iostream>
#include <iomanip>
#include <e_float/e_float.h>

int main(void)
{
  // Query the number of base-10 digits and print its value.
  static const int d = std::numeric_limits<e_float>::digits10;
  std::cout << d << std::endl;

  // Query the maximum value.
  static const e_float m = std::numeric_limits<e_float>::max();

  // Set the output stream precision and print the maximum value.
  std::cout << std::setprecision(std::numeric_limits<e_float>::digits10)
            << m << std::endl;
}

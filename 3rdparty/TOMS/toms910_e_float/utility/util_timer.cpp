#include <time.h>

#include <utility/util_timer.h>

double Util::timer::get_sec(void)
{
  static const double scale = static_cast<double>(CLOCKS_PER_SEC);
  const double tick = static_cast<double>(::clock());
  return static_cast<double>(tick / scale);
}

double Util::timer::elapsed(void) const
{
  const double delta = static_cast<double>(get_sec() - start);

  return static_cast<double>(delta - offset);
}

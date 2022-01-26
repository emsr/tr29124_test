#ifndef ROUNDING_TOOLS_H
#define ROUNDING_TOOLS_H 1

#include <cfenv>

/// The four different ronding modes.
enum struct fpround
{
  to_nearest  = FE_TONEAREST,
  upward      = FE_UPWARD,
  downward    = FE_DOWNWARD,
  toward_zero = FE_TOWARDZERO
};

void fpnear();
void fpdown();
void fpup();
void fpzero();

class round_sentinel
{
private:

  fpround m_round;

public:

  round_sentinel()
  { fpnear(); }

  round_sentinel(fpround rm)
  : m_round(static_cast<fpround>(std::fegetround()))
  {
    switch(rm)
    {
    case fpround::to_nearest:
      fpnear();
      break;
    case fpround::downward:
      fpdown();
      break;
    case fpround::upward:
      fpup();
      break;
    case fpround::toward_zero:
      fpzero();
      break;
    }
  }

  ~round_sentinel()
  { fesetround(static_cast<int>(this->m_round)); }
};

inline void
fpnear()
{
  fesetround(FE_TONEAREST);
}

inline void
fpdown()
{
  fesetround(FE_DOWNWARD);
}

inline void
fpup()
{
  fesetround(FE_UPWARD);
}

inline void
fpzero()
{
  fesetround(FE_TOWARDZERO);
}

inline float
to_float(const double& d, fpround rm)
{
  round_sentinel rs(rm);

  auto fres = float(d);

  return fres;
}

#endif // ROUNDING_TOOLS_H

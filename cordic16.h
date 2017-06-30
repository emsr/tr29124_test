// Cordic in 16-bit signed fixed point math
// Function is valid for arguments in range -pi/2 -- pi/2
// for values pi/2 -- pi: value = cordic16_half_pi - (theta - cordic16_half_pi)
// and similarly for values -pi -- -pi/2
//
// 1.0 = 16384
// 1/K = 0.6072529350088812561694
// pi = 3.1415926536897932384626
//Constants
constexpr short cordic16_1K = 0x26DD;
constexpr short cordic16_half_pi = 0x6487;
//#define MUL 16384.0
constexpr short cordic16_ntab = 16;
constexpr short
cordic16_ctab[cordic16_ntab]
{0x3243, 0x1DAC, 0x0FAD, 0x07F5,
 0x03FE, 0x01FF, 0x00FF, 0x007F,
 0x003F, 0x001F, 0x000F, 0x0007,
 0x0003, 0x0001, 0x0000, 0x0000,};

std::pair<short, short>
cordic(short theta, int n)
{
  short x = cordic16_1K;
  short y = 0;
  short z = theta;
  n = (n > cordic16_ntab) ? cordic16_ntab : n;
  for (int k = 0; k < n; ++k)
    {
      // Get sign.
      short d = z >= 0 ? 0 : -1;
      short tx = x - (((y >> k) ^ d) - d);
      short ty = y + (((x >> k) ^ d) - d);
      z -= ((cordic16_ctab[k] ^ d) - d);
      x = tx;
      y = ty;
    }  
  return std::make_pair(y, x);
}


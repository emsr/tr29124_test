// Cordic in 32-bit signed fixed point math
// Function is valid for arguments in range -pi/2 -- pi/2
// for values pi/2 -- pi: value = cordic32_half_pi - (theta-cordic32_half_pi)
// and similarly for values -pi -- -pi/2
//
// 1.0 = 1073741824 // 2147483647 >> 2 + 1?
// 1/K = 0.6072529350088812561694
// pi = 3.1415926535897932384626
//Constants
constexpr int cordic32_1K = 0x26DD3B6A;
constexpr int cordic32_half_pi = 0x6487ED51;
constexpr int cordic32_one = 1073741824;
constexpr int cordic32_ntab = 32;
constexpr int
cordic32_atan_tab[cordic32_ntab]
{
  0x3243F6A8, 0x1DAC6705, 0x0FADBAFC, 0x07F56EA6,
  0x03FEAB76, 0x01FFD55B, 0x00FFFAAA, 0x007FFF55,
  0x003FFFEA, 0x001FFFFD, 0x000FFFFF, 0x0007FFFF,
  0x0003FFFF, 0x0001FFFF, 0x0000FFFF, 0x00007FFF,
  0x00003FFF, 0x00001FFF, 0x00000FFF, 0x000007FF, 
  0x000003FF, 0x000001FF, 0x000000FF, 0x0000007F,
  0x0000003F, 0x0000001F, 0x0000000F, 0x00000008,
  0x00000004, 0x00000002, 0x00000001, 0x00000000
};

short
atan2_cordic(int n, short y0, short x0)
{
  short x = x0;
  short y = y0;
  short z = 0;
  for (int k = 0; k < n; ++k)
    {
      // Get sign.
      short d = y > 0 ? -1 : 0;
      short tx = x - (((y >> k) ^ d) - d);
      short ty = y + (((x >> k) ^ d) - d);
      z -= ((cordic32_atan_tab[k] ^ d) - d);
      x = tx;
      y = ty;
    }
  return z;
}

std::pair<int, int>
sincos_cordic(int n, int theta)
{
  int x = cordic32_1K;
  int y = 0;
  int z = theta;
  n = n > cordic32_ntab ? cordic32_ntab : n;
  for (int k = 0; k < n; ++k)
    {
      // Get sign.
      int d = z >= 0 ? 0 : -1;
      int tx = x - (((y >> k) ^ d) - d);
      int ty = y + (((x >> k) ^ d) - d);
      z -= ((cordic32_atan_tab[k] ^ d) - d);
      x = tx;
      y = ty;
    }
  return std::make_pair(y, x);
}

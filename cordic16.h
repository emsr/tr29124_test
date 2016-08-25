//Cordic in 16 bit signed fixed point math
//Function is valid for arguments in range -pi/2 -- pi/2
//for values pi/2--pi: value = half_pi-(theta-half_pi) and similarly for values -pi---pi/2
//
// 1.0 = 16384
// 1/k = 0.6072529350088812561694
// pi = 3.1415926536897932384626
//Constants
#define cordic_1K 0x26DDs
#define half_pi 0x6487s
#define MUL 16384.0
#define CORDIC_NTAB 16
short
cordic_ctab[16]{0x3243s, 0x1DACs, 0x0FADs, 0x07F5s,
                0x03FEs, 0x01FFs, 0x00FFs, 0x007Fs,
	        0x003Fs, 0x001Fs, 0x000Fs, 0x0007s,
	        0x0003s, 0x0001s, 0x0000s, 0x0000s,};

void
cordic(short theta, short *s, short *c, int n)
{
  short x = cordic_1K;
  short y = 0;
  short z = theta;
  int n = (n > CORDIC_NTAB) ? CORDIC_NTAB : n;
  for (int k = 0; k < n; ++k)
    {
      // Get sign.
      short d = z >= 0 ? 0 : -1;
      short tx = x - (((y >> k) ^ d) - d);
      short ty = y + (((x >> k) ^ d) - d);
      short tz = z - ((cordic_ctab[k] ^ d) - d);
      x = tx;
      y = ty;
      z = tz;
    }  
  *c = x;
  *s = y;
}


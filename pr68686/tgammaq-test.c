// Link with -lm -lquadmath

#include <stdio.h>
#include <math.h>
#include <quadmath.h>

int main(void)
{
	FILE *f1, *f2;
        int n;

        f1 = fopen("tgamma.dat", "w");
        if (f1) {
                for (n = -500; n <= 500; ++n) {
                        double x = n / 100.0;
                        double y = tgamma(x);
                        if (!isnan(y)) {
                                fprintf(f1, "%.16f %.16f\n", x, y);
                        }
                }
                fclose(f1);
        }

        f2 = fopen("tgammaq.dat", "w");
        if (f2) {
                for (n = -500; n <= 500; ++n) {
                        __float128 x = n / 100.0Q;
                        __float128 y = tgammaq(x);
                        if (!isnanq(y)) {
                                char buf[256];
                                quadmath_snprintf(buf, sizeof buf, "%.33Qf", x);
                                fputs(buf, f2);
				fputc(' ', f2);
                                quadmath_snprintf(buf, sizeof buf, "%.33Qf", y);
                                fputs(buf, f2);
				fputc('\n', f2);
                        }
                }
                fclose(f2);
        }

        return 0;
}

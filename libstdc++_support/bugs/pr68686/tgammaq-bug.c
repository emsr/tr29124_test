#include <stdio.h>
#include <quadmath.h>

int main(void)
{
	// Expected value: 2.36327...
	// Actual value: -2.36327...
	__float128 result = tgammaq(-1.5Q);
	char buf[256];
	quadmath_snprintf(buf, sizeof buf, "%.33Qf", result);
	puts(buf);

	result = tgammaq(-2.5Q);
	quadmath_snprintf(buf, sizeof buf, "%.33Qf", result);
	puts(buf);

	return 0;
}

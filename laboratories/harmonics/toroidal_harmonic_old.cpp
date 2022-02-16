/*
g++ -Wall -Wextra -o tor_old toroidal_harmonic_old.cpp
./tor_old > tor_old.txt
*/

#include <cmath>
#include <iostream>

double gammah(int n, double huge);
double frac(double z, int m, int n_max, double eps, double sqrttiny);
double fracps(double qz, int m, int n, double eps, double sqrttiny);
double factco(int n, double pl, int m);
double expan(double z, int mode, int iprec, double huge, double qargu, int m);
double psi(int m, int k);
double factor(int m, int k);

/**
 *  INPUT :
 * @param z        Argument of the functions
 *    DimM  M-dimension of the arrays. DimM must be greater than m_max.
 *    DimN  N-dimension of the arrays. DimN must be greater than n_max.
 * @param m_max Maximum order of the functions.
 *          We calculate functions of all orders below m_max.
 * @param n_max Maximum degree of the functions.
 *          We calculate functions of all degrees below min(n_max, n_new).
 *
 *   *If mode IS EQUAL TO 0:
 *    Pl[m][n]
 *             THESE VALUES ARE KEPT IN AN ARRAY
 *    Ql[m][n]
 *             THESE VALUES ARE KEPT IN AN ARRAY
 *
 *    m_new     MAXIMUM  ORDER OF FUNCTIONS CALCULATED WHEN
 *             Ql (m_max,0)   IS LARGER THAN 1/tiny
 *              (OVERFLOW LIMIT = 1/tiny, tiny IS DEFINED BELOW).
 *    n_new     MAXIMUM  DEGREE OF FUNCTIONS CALCULATED WHEN
 *             Pl (m,n_max)   IS LARGER THAN 1/tiny  FOR SOME
 *             m=0,...,m_new
 *              (OVERFLOW LIMIT = 1/tiny, tiny IS DEFINED BELOW).
 *    NOTE1: FOR A PRECISION OF 10**(-12), if z>5 AND (z/m)>0.22
 *           THE CODE USES A SERIES EXPANSION FOR Pl[m][0]. WHEN
 *           z<20 AND (z/m)<0.22 A CONTINUED FRACTION IS APPLIED.
 *    NOTE2: FOR A PRECISION OF 10**(-8), if z>5 AND (z/m)>0.12
 *           THE CODE USES A SERIES EXPANSION FOR Pl[m][0]. WHEN
 *           z<20 AND (z/m)<0.12 A CONTINUED FRACTION IS APPLIED.
 *
 *   *If mode IS EQUAL TO 1:
 *      THE SET OF FUNCTIONS EVALUATED IS:
 *                 Pl[m][n]/gamma(m+1/2),Ql[m][n]/gamma(m+1/2),
 *      WHICH ARE RESPECTIVELY STORED IN THE ARRAYS Pl[m][n],Ql[m][n]
 *      m_new AND n_new REFER TO THIS NEW SET OF FUNCTIONS
 *      NOTE1 AND NOTE2 ALSO APPLY IN THIS CASE
 *   *If mode IS EQUAL TO 2:
 *      THE CODE PERFORMS AS FOR mode 1, BUT THE RESTRICTION z<20
 *      FOR THE EVALUATION OF THE CONTINUED FRACTION IS NOT
 *      CONSIDERED
 *      WARNING: USE ONLY if HIGH m'S FOR z>20 ARE REQUIRED. THE
 *      EVALUATION OF THE CF MAY FAIL TO CONVERGE FOR TOO HIGH z'S
 *  PARAMETERS:
 *   mode: ENABLES THE ENLARGEMENT OF THE RANGE OF ORDERS AND
 *         DEGREES THAT CAN BE EVALUATED.
 *   s_eps:  CONTROLS THE ACCURACY OF THE CONTINUED FRACTIONS
 *         AND SERIES.
 *   iprec: REQUIRED PRECISION IN THE EVALUATION OF TOROIDAL
 *         HARMONICS.
 *           *If iprec=1, PRECISION=10**(-12) (TAKING s_eps<10**(-12))
 *           *If iprec=2, PRECISION=10**(-8) (TAKING s_eps<10**(-8))
 *   tiny: SMALL PARAMETER NEAR THE UNDERFLOW LIMIT.
 */
template <int DimM, int DimN>
void
DTORH3(double z, const int m_max, const int n_max, double Pl[DimM][DimN], double Ql[DimM][DimN], int &m_new, int &n_new)
{
    int mp, ical, m, n, np, i, l;
    double qz, qdc1, qargu, argu1, d1, d2, fl, cc, ar, gamma, dfacqs, fcp, dd, pl0, qm0, dfac3, fc, dfacc, fc2, dfac4, ELLIP1,
        ELLIP2, pr[3];
    constexpr double pi = 3.14159265358979323, s_eps = 1.0e-14, tiny = 1.0e-290;
    constexpr int mode = 0, iprec = 1;
    constexpr auto huge = 1.0 / tiny;
    const auto sqrttiny = std::sqrt(tiny);
    if (iprec != 1 && iprec != 2)
    {
        std::cout << "iprec must be 1 or 2";
        return;
    }
    pr[0] = 0.0; // Dummy - was 1-offset.
    pr[1] = 0.22;
    pr[2] = 0.12;
    // s_eps: Required accuracy for the continued fraction (modified Lentz)
    // tiny: Small parameter to prevent overflows in the CF (close to the underflow limit)
    if (z <= 1.0)
    {
        std::cout << "Improper argument: z must be greater than 1";
        return;
    }
    qz = z;
    qdc1 = qz * qz - 1.0;
    qargu = qz / std::sqrt(qdc1);
    argu1 = std::sqrt(2.0 / (z + 1.0));
    const auto sqrtpi = std::sqrt(pi);

    // We evaluate Q^{(0)}_{-1/2},Q^{(1)}_{-1/2}
    // Using SLATEC routines for elliptic functions
    Ql[0][0] = argu1 * std::comp_ellint_1(argu1);
    Ql[1][0] = -1.0 / std::sqrt(2.0 * (qz - 1.0)) * std::comp_ellint_2(argu1);

    // we apply forward recurrence in m for Q's
    mp = 1;
    if (mode == 0)
    {
        while ((mp <= m_max) && (std::abs(Ql[mp][0]) < huge))
        {
            Ql[mp + 1][0] = -2.0 * mp * qargu * Ql[mp][0] - (mp - 0.5) * (mp - 0.5) * Ql[mp - 1][0];
            mp = mp + 1;
        }
        if ((mp - 1) < m_max)
            m_max = mp - 1;
        m_new = m_max;
    }
    else
    {
        Ql[0][0] = Ql[0][0] / sqrtpi;
        Ql[1][0] = Ql[1][0] * 2.0 / sqrtpi;
        while ((mp <= m_max) && (std::abs(Ql[mp][0]) < huge))
        {
            d1 = mp + 0.5;
            Ql[mp + 1][0] = -2.0 * mp * qargu * Ql[mp][0] / d1 - (mp - 0.5) * Ql[mp - 1][0] / d1;
            mp = mp + 1;
        }
        if ((mp - 1) < m_max)
            m_max = mp - 1;
        m_new = m_max;
    }
    fl = m_max / 2.0;
    cc = std::abs(double(int(fl)) - fl);
    if (cc < 0.4)
    {
        ar = 1.0;
    }
    else
    {
        ar = -1.0;
    }
    if (mode == 0)
    {
        gamma = gammah(m_max, huge) * ar * sqrtpi;
        if (std::abs(gamma) < tiny)
        {
            std::cout << "m_max is too large for mode=0";
            std::cout << "Better try mode=1";
            return;
        }
    }
    else
    {
        gamma = ar;
    }
    dfacqs = -(gamma / Ql[m_max + 1][0]) * gamma / pi;
    // EVALUATION OF Pl[m_max][0],Pl[m_max+1][0]
    // WE CHOOSE EXPANSION OR CF FOR Pl[m_max][0]
    // DEPENDING ON THE VALUES OF z,m_max AND mode
    ical = 1;
    if ((z / m_max) > pr[iprec])
        ical = 2;
    if (z < 5.0)
    {
        ical = 1;
    }
    else if (z > 20.0)
    {
        if ((mode != 2) && (ical == 1))
            ical = 0;
    }
    if (ical == 0)
    {
        std::cout << "You must choose mode=2";
        return;
    }
    if (ical == 1)
    {
        // We calculate the CF for P's
        fcp = fracps(qargu, m_max + 1, 0, s_eps, sqrttiny);
        dd = m_max + 0.5;
        if (mode != 0)
        {
            fcp = fcp / dd;
            dfacqs = dfacqs / dd;
        }
        Pl[m_max][0] = dfacqs / std::sqrt(qdc1) / (1.0 - fcp * Ql[m_max][0] / Ql[m_max + 1][0]);
        Pl[m_max + 1][0] = fcp * Pl[m_max][0];
    }
    else
    {
        pl0 = expan(z, mode, iprec, huge, qargu, m_max);
        Pl[m_max][0] = pl0;
        dd = m_max + 0.5;
        if (mode != 0)
        {
            fcp /= dd;
            dfacqs /= dd;
        }
        Pl[m_max + 1][0] = (Ql[m_max + 1][0] / Ql[m_max][0]) * (Pl[m_max][0] - dfacqs / std::sqrt(qdc1));
    }
    // Evaluation of Pl[m_max][1]
    qm0 = Ql[m_max][0];
    dfac3 = (gamma / qm0) * gamma / pi / (0.5 - m_max);
    fc = frac(z, m_max, 0, s_eps, sqrttiny);
    Pl[m_max][1] = Pl[m_max][0] * fc + dfac3;
    // EVALUATION OF Pl[m_max+1][1]
    qm0 = Ql[m_max + 1][0];
    if (mode == 0)
    {
        dfacc = (0.5 + m_max);
    }
    else
    {
        dfacc = 1.0 / (0.5 + m_max);
    }
    dfac3 = -(gamma / qm0) * gamma * dfacc / pi;
    fc2 = frac(z, m_max + 1, 0, s_eps, sqrttiny);
    Pl[m_max + 1][1] = Pl[m_max + 1][0] * fc2 + dfac3;
    // We apply backward recurrence over m to get the set Pl[m][0],Pl[m][1]
    if (mode == 0)
    {
        for (int i = 1; i <= m_max; ++i)
        {
            mp = m_max + 1 - i;
            m = mp - 1;
            Pl[m][0] = -(Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / ((0.5 - mp) * (0.5 - mp));
        }
        for (int i = 1; i <= m_max; ++i)
        {
            mp = m_max + 1 - i;
            m = mp - 1;
            Pl[m][1] = (Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) / ((1.5 - mp) * (0.5 + mp));
        }
    }
    else
    {
        for (int i = 1; i <= m_max; ++i)
        {
            mp = m_max + 1 - i;
            m = mp - 1;
            Pl[m][0] = ((mp + 0.5) * Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / (0.5 - mp);
        }
        for (int i = 1; i <= m_max; ++i)
        {
            mp = m_max + 1 - i;
            m = mp - 1;
            Pl[m][1] = ((mp + 0.5) * Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) * (mp - 0.5) / ((1.5 - mp) * (0.5 + mp));
        }
    }
    // Now, we perform the evaluation over n for each m
    // We start evaluating the set of Pl[m][n]
    for (int l = 0; l <= m_max; ++l)
    {
        m = l;
        fl = m / 2.0;
        cc = std::abs(double(int(fl)) - fl);
        if (cc < 0.4)
        {
            ar = 1.0;
        }
        else
        {
            ar = -1.0;
        }
        if (mode == 0)
        {
            gamma = gammah(m, huge) * ar * sqrtpi;
        }
        else
        {
            gamma = ar;
        }
        np = 1;
        while ((np <= n_max) && (std::abs((np - m + 0.5) * Pl[l][np]) < huge))
        {
            Pl[l][np + 1] = (2.0 * np * z * Pl[l][np] - (np + m - 0.5) * Pl[l][np - 1]) / (np - m + 0.5);
            np = np + 1;
        }
        n_max = np - 1;
        n_new = n_max;
    }
    // We evaluate the set of Ql[m][n]
    for (int l = 0; l <= m_max; ++l)
    {
        m = l;
        fl = m / 2.0;
        cc = std::abs(double(int(fl)) - fl);
        if (cc < 0.4)
        {
            ar = 1.0;
        }
        else
        {
            ar = -1.0;
        }
        if (mode == 0)
        {
            gamma = gammah(m, huge) * ar * sqrtpi;
        }
        else
        {
            gamma = ar;
        }
        // We evaluate the c.f. for Q's using Lentz-Thompson
        // If m=0,1. In other case we calculate Ql[m][n_max] and
        // Ql[m][n_max+1] from the recurrence.
        if (m == 2)
        {
            if (mode == 0)
            {
                mp = 1;
                while (mp <= m_max)
                {
                    Ql[mp + 1][n_max]
                        = -2.0 * mp * qargu * Ql[mp][n_max] - (mp - n_max - 0.5) * (mp + n_max - 0.5) * Ql[mp - 1][n_max];
                    mp = mp + 1;
                }
                mp = 1;
                while (mp <= m_max)
                {
                    Ql[mp + 1][n_max + 1] = -2.0 * mp * qargu * Ql[mp][n_max + 1]
                                            - (mp - n_max - 1.5) * (mp + n_max + 0.5) * Ql[mp - 1][n_max + 1];
                    mp = mp + 1;
                }
            }
            else
            {
                mp = 1;
                while (mp <= m_max)
                {
                    d1 = mp + 0.5;
                    d2 = mp - 0.5;
                    Ql[mp + 1][n_max] = -2.0 * mp * qargu * Ql[mp][n_max] / d1
                                        - (mp - n_max - 0.5) * (mp + n_max - 0.5) * Ql[mp - 1][n_max] / (d1 * d2);
                    mp = mp + 1;
                }
                mp = 1;
                while (mp <= m_max)
                {
                    d1 = mp + 0.5;
                    d2 = mp - 0.5;
                    Ql[mp + 1][n_max + 1] = -2.0 * mp * qargu * Ql[mp][n_max + 1] / d1
                                            - (mp - n_max - 1.5) * (mp + n_max + 0.5) * Ql[mp - 1][n_max + 1] / (d1 * d2);
                    mp = mp + 1;
                }
            }
        }
        if ((m == 0) || (m == 1))
        {
            fc = frac(z, m, n_max, s_eps, sqrttiny);
            dfac4 = factco(n_max, Pl[l][n_max + 1], m) * gamma * gamma / pi / (n_max + m + 0.5);
            // Evaluation of Ql[l][n_max+1] AND Ql[l][n_max] USING
            // The Wronskian W{Pl[l][n_max],Ql[l][n_max]},
            // The known values of Pl[l][n_max+1] AND Pl[l][n_max]
            // The value of H = Ql[l][n_max+1]/Ql[l][n_max]
            Ql[l][n_max] = dfac4 / (1.0 - fc * Pl[l][n_max] / Pl[l][n_max + 1]);
            Ql[l][n_max + 1] = Ql[l][n_max] * fc;
        }
        // We use the backward recurrence relation for Q's
        for (int i = 1; i <= n_max; ++i)
        {
            np = n_max + 1 - i;
            n = np - 1;
            Ql[l][n] = ((np + np) * z * Ql[l][np] - (np - m + 0.5) * Ql[l][np + 1]) / (np + m - 0.5);
        }
    }
    return;
}

double
gammah(int n, double over)
{
    int i, j;
    i = n;
    j = 2 * i - 1;
    double gh = 1.0;
    while ((j >= 1) && (gh < over))
    {
        gh *= j / 2.0;
        i = i - 1;
        j = 2 * i - 1;
    }
    if (j > 1)
    {
        gh = 0.0;
    }
    return gh;
}

double
frac(double z, int m, int n_max, double eps, double sqrttiny)
{
    const int n = n_max;
    const double dz2 = z + z;
    const double dn0 = n + m;
    const double dn1 = n + 1.0;
    const double dn2 = dn0 + 0.5;
    const double dn3 = dn0 - 0.5;
    const double dn4 = 1.0 - m - m;
    double b = 2.0 * dn1 * z / dn2;
    double a = 1.0;
    double fc = sqrttiny;
    double c0 = fc;
    double d0 = 0.0;

    int mm = 0;
    while (mm < 10000)
    {
        d0 = b + a * d0;
        if (std::abs(d0) < sqrttiny)
            d0 = sqrttiny;
        c0 = b + a / c0;
        if (std::abs(c0) < sqrttiny)
            c0 = sqrttiny;
        d0 = 1.0 / d0;
        double delta = c0 * d0;
        fc *= delta;
        ++mm;
        a = -(1.0 + dn4 / (dn3 + mm));
        b = dz2 * (dn1 + mm) / (dn2 + mm);
        if (std::abs(delta - 1.0) < eps)
            break;
    }
    if (mm == 10000)
    {
        std::cout << "CF convergence fails\n";
        return fc;
    }
    return fc;
}

double
fracps(double qz, int m, int n, double eps, double sqrttiny)
{
    const double dn2 = n * n;
    const double dz2 = qz + qz;
    const double dm = m - 0.5;
    double b = dz2 * m;
    double a = dn2 - dm * dm;
    double fc = sqrttiny;
    double c0 = fc;
    double d0 = 0.0;

    int mm = 0;
    while (mm < 10000)
    {
        d0 = b + a * d0;
        if (std::abs(d0) < sqrttiny)
            d0 = sqrttiny;
        c0 = b + a / c0;
        if (std::abs(c0) < sqrttiny)
            c0 = sqrttiny;
        d0 = 1.0 / d0;
        double delta = c0 * d0;
        fc *= delta;
        ++mm;
        a = dn2 - (mm + dm) * (mm + dm);
        b = dz2 * (mm + m);
        if (std::abs(delta - 1.0) < eps)
            break;
    }
    if (mm == 10000)
    {
        std::cout << "CF convergence fails\n";
        return fc;
    }
    return fc;
}

double
factco(int n, double pl, int m)
{
    double fact = 1.0 / pl;
    double x1 = m + 0.5;
    double x2 = -m + 0.5;
    int j = n;
    while (j >= 0)
    {
        fact *= (j + x1) / (j + x2);
        --j;
    }
    return fact;
}

double
expan(double z, int mode, int iprec, double huge, double qargu, int m)
{
    constexpr double pi = 3.14159265358979323;
    double preci[3];
    preci[0] = 0.0; // Dummy.
    preci[1] = 1.0e-13;
    preci[2] = 1.0e-9;
    const double sqrtpi = std::sqrt(pi);
    const double db = 2.0 * std::log(2.0);
    const double fl = m / 2.0;
    const double cc = std::abs(double(int(fl)) - fl);
    double ar;
    if (cc < 0.4)
        ar = 1.0;
    else
        ar = -1.0;
    double dz = 1.0;
    for (int i = 1; i <= m; ++i)
        dz /= qargu;
    double gamma;
    if (mode == 0)
        gamma = gammah(m, huge) * ar * sqrtpi;
    else
        gamma = ar;
    const double dfac = 2.0 / pi * dz * gamma / sqrtpi;
    const double df1 = std::log(2.0 * z);
    double a0 = 1.0 / std::sqrt(2.0 * z);
    const double z2i = 1.0 / (z * z);
    double delta = 1.0;
    double sum = 0.0;
    double da2 = factor(m, 0);
    double da1 = db + psi(m, 0);

    int k = 0;
    while (k < 1000)
    {
        delta = (df1 + da1) * da2 * a0;
        sum += delta;
        double dcc = 0.5 + m + 2.0 * k;
        double dccp = dcc + 1.0;
        double dkk = k + 1.0;
        da2 = da2 * dccp * dcc / (dkk * dkk);
        da1 = da1 + 1.0 / dkk - 1.0 / dccp - 1.0 / dcc;
        ++k;
        a0 *= 0.25 * z2i;
        if (std::abs(delta / sum) < preci[iprec])
            break;
    }
    double pl = sum * dfac;
    return pl;
}

double
psi(int m, int k)
{
    double factr1 = 0.0;
    double factr2 = 0.0;
    int n = 2 * k + m;

    int j = 1;
    if (k >= 1)
    {
        while (j <= k)
        {
            factr1 += 1.0 / j;
            ++j;
        }
    }

    int i = 1;
    while (i <= n)
    {
        factr2 += 1.0 / (2.0 * i - 1.0);
        ++i;
    }
    if (n == 0)
        factr2 = 0.0;

    return factr1 - 2.0 * factr2;
}

double
factor(int m, int k)
{
    double fact = 1.0;
    if (k >= 1)
    {
        const double x1 = m + 0.5;
        const int n = 2 * k - 1;
        int i = k;
        int j = n;
        while (j >= 0)
        {
            fact *= (j + x1) / i / i;
            j = j - 1;
            i = i - 1;
            if (i == 0)
                i = 1;
        }
    }
    return fact;
}

int
main()
{
}

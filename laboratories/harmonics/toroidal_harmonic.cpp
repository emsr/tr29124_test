/*
g++ -Wall -Wextra -o toroidal_harmonic toroidal_harmonic.cpp
./toroidal_harmonic > toroidal_harmonic.txt
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

double gammah(int n, double huge);
double frac(double z, int m, int n_max, double eps, double sqrttiny);
double fracps(double qz, int m, int n, double eps, double sqrttiny);
double factco(int n, double pl, int m);
double expan(double z, int mode, int iprec, double huge, double qargu, int m);
double psi(int m, int k);
double factor(int m, int k);
void toroidal_asymp_series(double a, int m, double &sa, double &sb);
double toroidal_asymp(double z, int m, int mode, double gammapr);
void
toroidal_harmonic_dual(double z, int m_max, int n_max, std::vector<std::vector<double>> &Pl,
                       std::vector<std::vector<double>> &Ql, int &m_new, int &n_new);
void
toroidal_harmonic_primal(double z, int m_max, int n_max, std::vector<std::vector<double>> &Pl,
                         std::vector<std::vector<double>> &Ql, int &m_new, int &n_new);

void
toroidal_harmonic(double z, int m_max, int n_max, std::vector<std::vector<double>> &Pl,
                  std::vector<std::vector<double>> &Ql, int &m_new, int &n_new)
{
    constexpr double sqrt2 = 1.4142135623730950488;

    if (z <= 1.0)
    {
        std::cout << "Improper argument: z must be greater than 1";
        return;
    }

    if (z < sqrt2)
      toroidal_harmonic_dual(z, m_max, n_max, Pl, Ql, m_new, n_new);
    else
      toroidal_harmonic_primal(z, m_max, n_max, Pl, Ql, m_new, n_new);
}

void
toroidal_harmonic_dual(double z, int m_max, int n_max, std::vector<std::vector<double>> &Pl,
                       std::vector<std::vector<double>> &Ql, int &m_new, int &n_new)
{
    constexpr double pi = 3.1415926535897932385;
    constexpr double sqrt2 = 1.4142135623730950488;
    // Controls the accuracy of the continued fractions and series.
    constexpr double s_eps = 1.0e-14;
    // Small parameter near the underflow limit.
    constexpr double tiny = 1.0e-290;
    // Enables the enlargement of the range of orders and degrees that can be evaluated.
    constexpr int mode = 0;
    // Required precision in the evaluation of toroidal harmonics.
    // If iprec=1, precision=10^{-12} (taking s_eps < 10^{-12})
    // If iprec=2, precision=10^{-8} (taking s_eps < 10^{-8})
    constexpr int iprec = 1;
    constexpr auto huge = 1.0 / tiny;
    constexpr auto hugeq = 1.0e-30 * huge;
    const auto sqrttiny = std::sqrt(tiny);
    double pr[3];
    pr[0] = 0.0; // Dummy - was 1-offset.
    pr[1] = 0.22;
    pr[2] = 0.12;

    // Dual algorithm
    // We evaluate P^{0}_{-1/2},P^{0}_{+1/2}
    const auto qz = z;
    const auto qdc1 = qz * qz - 1.0;
    const auto qq = std::sqrt(qdc1);
    const auto qargu = qz / qq;
    const auto xl = qargu;
    const auto dzz = std::sqrt(z * z - 1.0);
    const auto argu1 = std::sqrt((z - 1.0) / (z + 1.0));
    const auto argu2 = std::sqrt(2.0 * dzz / (z + dzz));
    const auto sqrtpi = std::sqrt(pi);
    const auto pi3d2 = sqrtpi * pi;
    Pl[0][0] = 2.0 / pi * std::sqrt(2.0 / (z + 1.0)) * std::comp_ellint_1(argu1);
    Pl[0][1] = 2.0 / pi * std::sqrt(z + dzz) * std::comp_ellint_2(argu2);
    // We apply forward recurrence in n for P's
    if (mode == 0)
    {
        int np = 1;
        while (np <= n_max && abs(Pl[0][np]) < huge)
        {
            Pl[0][np + 1] = (2.0 * np * z * Pl[0][np] - (np - 0.5) * Pl[0][np - 1]) / (np + 0.5);
            ++np;
        }
        if (np - 1 != n_max)
            n_max = np - 1;
        n_new = n_max;
    }
    else
    {
        Pl[0][0] /= sqrtpi;
        Pl[0][1] /= sqrtpi;
        int np = 1;
        while (np <= n_max && abs(Pl[0][np]) < huge)
        {
            Pl[0][np + 1] = (2.0 * np * z * Pl[0][np] - (np - 0.5) * Pl[0][np - 1]) / (np + 0.5);
            ++np;
        }
        if ((np - 1) != n_max)
            n_max = np - 1;
        n_new = n_max;
    }
    // We choose cf, series or asymptotic expansion
    // depending on the values of xl and n_max
    int ical = 1;
    if (xl / n_max > pr[iprec])
        ical = 2;
    if (xl < 5.0)
        ical = 1;
    else if (xl > 20.0 && ical == 1)
        ical = 0;

    if (ical == 1)
    {
        // We evaluate the continued fraction
        const auto fc = frac(z, 0, n_max, s_eps, sqrttiny);
        // We evaluate Ql(0,n_max+1) and Ql(0,n_max) using :
        // the Wronskian : W{Pl(0,n_max),Ql(0,n_max)} =1./(1.-z**2)
        // the known values of Pl(0,n_max+1) AND Pl(0,n_max)
        // the value of fc=Ql(0,n_max+1)/Ql(0,n_max)
        if (mode == 0)
            Ql[0][n_max] = 1.0 / (Pl[0][n_max + 1] - fc * Pl[0][n_max]) / (n_max + 0.5);
        else
            Ql[0][n_max] = 1.0 / (Pl[0][n_max + 1] - fc * Pl[0][n_max]) / ((n_max + 0.5) * pi);
        Ql[0][n_max + 1] = Ql[0][n_max] * fc;
    }
    else if (ical == 2)
    {
        // We calculate series
        _40:
        const auto fl = n_max / 2.0;
        const auto cc = std::abs(double(int(fl)) - fl);
        double ar;
        if (cc < 0.4)
            ar = 1.0;
        else
            ar = -1.0;

        const auto gammapr = ar * gammah(n_max, huge);
        const auto gamma = gammapr * sqrtpi;
        if (std::abs(gamma) < tiny)
        {
            n_max -= 5;
            goto _40;
        }
        const auto dfac = pi3d2 / sqrt2 * std::pow(xl * xl - 1.0, 0.25);
        const auto pl0 = expan(xl, 0, iprec, huge, z, n_max);
        Ql[0][n_max] = (pl0 / gamma) * dfac;
        if (mode == 0)
            Ql[0][n_max + 1] = (Pl[0][n_max + 1] * Ql[0][n_max] - 1.0 / (n_max + 0.5)) / Pl[0][n_max];
        else
        {
            Ql[0][n_max] = Ql[0][n_max] / sqrtpi;
            Ql[0][n_max + 1] = (Pl[0][n_max + 1] * Ql[0][n_max] - 1.0 / (n_max + 0.5) / pi) / Pl[0][n_max];
        }
    }
    else if (ical == 0)
    {
        // We use the asymptotic expansion
        _50:
        const auto fl = n_max / 2.0;
        const auto cc = std::abs(double(int(fl)) - fl);
        double ar;
        if (cc < 0.4)
            ar = 1.0;
        else
            ar = -1.0;
        const auto gammapr = ar * gammah(n_max, huge);
        const auto gamma = gammapr * sqrtpi;
        if (std::abs(gamma) < tiny)
        {
            n_max = n_max - 5;
            goto _50;
        }
        const auto dfac = pi3d2 / sqrt2 * std::pow(xl * xl - 1.0, 0.25);
        const auto pl0 = toroidal_asymp(xl, n_max, 0, gammapr);
        Ql[0][n_max] = (pl0 / gamma) * dfac;
        if (mode == 0)
        {
            Ql[0][n_max + 1] = (Pl[0][n_max + 1] * Ql[0][n_max] - 1.0 / (n_max + 0.5)) / Pl[0][n_max];
        }
        else
        {
            Ql[0][n_max] = Ql[0][n_max] / sqrtpi;
            Ql[0][n_max + 1] = (Pl[0][n_max + 1] * Ql[0][n_max] - 1.0 / (n_max + 0.5) / pi) / Pl[0][n_max];
        }
    }
    // CFS for Ps and calculation of Pl(1,n_max+1),Pl(1,n_max)
    const auto fcp1 = fracps(xl, 1, n_max + 1, s_eps, sqrttiny);
    Pl[1][n_max + 1] = fcp1 * Pl[0][n_max + 1];
    const auto fcp2 = fracps(xl, 1, n_max, s_eps, sqrttiny);
    Pl[1][n_max] = fcp2 * Pl[0][n_max];
    if (mode == 1)
    {
        Pl[1][n_max + 1] = Pl[1][n_max + 1] * 2.0;
        Pl[1][n_max] = Pl[1][n_max] * 2.0;
    }
    // Evaluation of Ql(1,n_max+1),QL(1,n_max)
    const auto dizz = -1.0 / dzz;
    if (mode == 0)
    {
        Ql[1][n_max + 1] = (dizz + Pl[1][n_max + 1] * Ql[0][n_max + 1]) / Pl[0][n_max + 1];
        Ql[1][n_max] = (dizz + Pl[1][n_max] * Ql[0][n_max]) / Pl[0][n_max];
    }
    else
    {
        const auto da = 2.0 / pi;
        Ql[1][n_max + 1] = (dizz * da + Pl[1][n_max + 1] * Ql[0][n_max + 1]) / Pl[0][n_max + 1];
        Ql[1][n_max] = (dizz * da + Pl[1][n_max] * Ql[0][n_max]) / Pl[0][n_max];
    }
    // Evaluation of the set {Ql(m,n)}
    // We apply forward recurrence in m for Q's
    if (mode == 0)
    {
        int mp = 1;
        while (mp <= m_max && std::abs(Ql[mp][n_max]) < hugeq)
        {
            Ql[mp + 1][n_max]
                = -2.0 * mp * qargu * Ql[mp][n_max] - (mp - n_max - 0.5) * (mp + n_max - 0.5) * Ql[mp - 1][n_max];
            ++mp;
        }
        m_max = mp - 2;
        mp = 1;
        while (mp <= m_max && std::abs(Ql[mp][n_max + 1]) < hugeq)
        {
            Ql[mp + 1][n_max + 1]
                = -2.0 * mp * qargu * Ql[mp][n_max + 1] - (mp - n_max - 1.5) * (mp + n_max + 0.5) * Ql[mp - 1][n_max + 1];
            ++mp;
        }
        if (mp - 1 < m_max)
            m_max = mp - 1;
    }
    else
    {
        int mp = 1;
        while (mp <= m_max && std::abs(Ql[mp][n_max]) < huge)
        {
            const auto d1 = mp + 0.5;
            const auto d2 = mp - 0.5;
            Ql[mp + 1][n_max] = -2.0 * mp * qargu * Ql[mp][n_max] / d1
                                - (mp - n_max - 0.5) * (mp + n_max - 0.5) * Ql[mp - 1][n_max] / (d1 * d2);
            ++mp;
        }
        m_max = mp - 2;
        mp = 1;
        while (mp <= m_max && std::abs(Ql[mp][n_max + 1]) < huge)
        {
            const auto d1 = mp + 0.5;
            const auto d2 = mp - 0.5;
            Ql[mp + 1][n_max + 1] = -2.0 * mp * qargu * Ql[mp][n_max + 1] / d1
                                    - (mp - n_max - 1.5) * (mp + n_max + 0.5) * Ql[mp - 1][n_max + 1] / (d1 * d2);
            ++mp;
        }
        if (mp - 1 < m_max)
            m_max = mp - 1;
    }
    m_new = m_max + 1;

    // Finally, for each m=0,...,m_max applying recurrence backwards, obtain the set {Ql(m,n)}
    for (int l = 0; l <= m_max + 1; ++l)
    {
        for (int i = 1; i <= n_max; ++i)
        {
            int np = n_max + 1 - i;
            int n = np - 1;
            Ql[l][n] = ((np + np) * z * Ql[l][np] - (np - l + 0.5) * Ql[l][np + 1]) / (np + l - 0.5);
        }
    }
    // Evaluation of the set {Pl(m,n)}
    // Calculation of the cf for Ps
    const auto fl = m_max / 2.0;
    const auto cc = std::abs(double(int(fl)) - fl);
    double ar;
    if (cc < 0.4)
        ar = 1.0;
    else
        ar = -1.0;

    double gamma;
    if (mode == 0)
        gamma = gammah(m_max, huge) * ar * sqrtpi;
    else
        gamma = ar;

    auto dfacqs = -(gamma / Ql[m_max + 1][0]) * gamma / pi;
    auto fcp = fracps(xl, m_max + 1, 0, s_eps, sqrttiny);
    const auto dd = m_max + 0.5;
    if (mode != 0)
    {
        fcp /= dd;
        dfacqs /= dd;
    }
    Pl[m_max][0] = dfacqs / qq / (1.0 - fcp * Ql[m_max][0] / Ql[m_max + 1][0]);
    Pl[m_max + 1][0] = fcp * Pl[m_max][0];
    // Evaluation of Pl(m_max,1)
    auto qm0 = Ql[m_max][0];
    auto dfac3 = (gamma / qm0) * gamma / pi / (0.5 - m_max);
    const auto fc = frac(z, m_max, 0, s_eps, sqrttiny);
    Pl[m_max][1] = Pl[m_max][0] * fc + dfac3;
    // Evaluation of Pl(m_max+1,1)
    qm0 = Ql[m_max + 1][0];
    double dfacc;
    if (mode == 0)
        dfacc = (0.5 + m_max);
    else
        dfacc = 1.0 / (0.5 + m_max);
    dfac3 = -(gamma / qm0) * gamma * dfacc / pi;
    const auto fc2 = frac(z, m_max + 1, 0, s_eps, sqrttiny);
    Pl[m_max + 1][1] = Pl[m_max + 1][0] * fc2 + dfac3;
    // We apply backward recurrence over m to get the set Pl(m,0),Pl(m,1)
    if (mode == 0)
    {
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][0] = -(Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / ((0.5 - mp) * (0.5 - mp));
        }
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][1] = (Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) / ((1.5 - mp) * (0.5 + mp));
        }
    }
    else
    {
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][0] = ((mp + 0.5) * Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / (0.5 - mp);
        }
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][1] = ((mp + 0.5) * Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) * (mp - 0.5) / ((1.5 - mp) * (0.5 + mp));
        }
    }
    // Now, we perform the evaluation over n for each m
    // We start evaluating the set of Pl(m,n)
    for (int l = 0; l <= m_max + 1; ++l)
    {
        int m = l;
        const auto fl = m / 2.0;
        const auto cc = std::abs(double(int(fl)) - fl);
        double ar;
        if (cc < 0.4)
            ar = 1.0;
        else
            ar = -1.0;

        if (mode == 0)
            gamma = gammah(m, huge) * ar * sqrtpi;
        else
            gamma = ar;

        int np = 1;
        while (np <= n_max && std::abs((np - m + 0.5) * Pl[l][np]) < huge)
        {
            Pl[l][np + 1] = (2.0 * np * z * Pl[l][np] - (np + m - 0.5) * Pl[l][np - 1]) / (np - m + 0.5);
            ++np;
        }
        n_max = np - 1;
        n_new = n_max;
    }
}

/**
 * @param z        Argument of the functions
 * @param m_max Maximum order of the functions.
 *              We calculate functions of all orders below m_max.
 * @param n_max Maximum degree of the functions.
 *              We calculate functions of all degrees below min(n_max, n_new).
 *
 *   *If mode IS EQUAL TO 0:
 *    Pl[m][n] Array of P
 *    Ql[m][n] Array of Q
 *
 * @param m_new  Maximum order of functions calculated when
 *               Ql[m_max][0] is larger than 1/tiny
 *              (overflow limit = 1/tiny, tiny is defined below).
 * @param n_new  Maximum degree of functions calculated when
 *               Pl[m][n_max] is larger than 1/tiny for some
 *               m = 0, ... , m_new
 *              (overflow limit = 1/tiny, tiny is defined below).
 * Note1: For a precision of 10**(-12), if z>5 and (z/m)>0.22
 *        the code uses a series expansion for Pl[m][0]. When
 *        z<20 and (z/m) < 0.22 a continued fraction is applied.
 * Note2: For a precision of 10**(-8), if z>5 and (z/m)>0.12
 *        the code uses a series expansion for Pl[m][0]. When
 *        z<20 and (z/m) < 0.12 a continued fraction is applied.
 *
 * If mode is equal to 1:
 *   The set of functions evaluated is:
 *     Pl[m][n]/gamma(m+1/2), Ql[m][n]/gamma(m+1/2),
 *     which are respectively stored in the arrays Pl[m][n],Ql[m][n]
 *     m_new and n_new refer to this new set of functions
 *     note1 and note2 also apply in this case
 * If mode is equal to 2:
 *     The code performs as for mode 1, but the restriction z<20
 *     for the evaluation of the continued fraction is not
 *     considered
 *     Warning: Use only if high m's for z>20 are required. The
 *     evaluation of the CF may fail to converge for too high z's
 *  PARAMETERS:
 *   mode: Enables the enlargement of the range of orders and degrees that can be evaluated.
 *   iprec: Required precision in the evaluation of toroidal harmonics.
 */
void
toroidal_harmonic_primal(double z, int m_max, int n_max, std::vector<std::vector<double>> &Pl,
                         std::vector<std::vector<double>> &Ql, int &m_new, int &n_new)
{
    constexpr double pi = 3.1415926535897932385;
    // Controls the accuracy of the continued fractions and series.
    constexpr double s_eps = 1.0e-14;
    // Small parameter near the underflow limit.
    constexpr double tiny = 1.0e-290;
    // Enables the enlargement of the range of orders and degrees that can be evaluated.
    constexpr int mode = 0;
    // Required precision in the evaluation of toroidal harmonics.
    // If iprec=1, precision=10^{-12} (taking s_eps < 10^{-12})
    // If iprec=2, precision=10^{-8} (taking s_eps < 10^{-8})
    constexpr int iprec = 1;
    constexpr auto huge = 1.0 / tiny;
    const auto sqrttiny = std::sqrt(tiny);
    if (iprec != 1 && iprec != 2)
    {
        std::cout << "iprec must be 1 or 2";
        return;
    }
    double pr[3];
    pr[0] = 0.0; // Dummy - was 1-offset.
    pr[1] = 0.22;
    pr[2] = 0.12;
    // s_eps: Required accuracy for the continued fraction (modified Lentz)
    // tiny: Small parameter to prevent overflows in the CF (close to the underflow limit)
    const auto qz = z;
    const auto qdc1 = qz * qz - 1.0;
    const auto qq = std::sqrt(qdc1);
    const auto qargu = qz / qq;
    const auto argu1 = std::sqrt(2.0 / (z + 1.0));
    const auto sqrtpi = std::sqrt(pi);

    // We evaluate Q^{(0)}_{-1/2},Q^{(1)}_{-1/2}
    Ql[0][0] = argu1 * std::comp_ellint_1(argu1);
    Ql[1][0] = -1.0 / std::sqrt(2.0 * (qz - 1.0)) * std::comp_ellint_2(argu1);

    // We apply forward recurrence in m for Q's
    int mp = 1;
    if (mode == 0)
    {
        while (mp <= m_max && std::abs(Ql[mp][0]) < huge)
        {
            Ql[mp + 1][0] = -2.0 * mp * qargu * Ql[mp][0] - (mp - 0.5) * (mp - 0.5) * Ql[mp - 1][0];
            ++mp;
        }
        if ((mp - 1) < m_max)
            m_max = mp - 1;
        m_new = m_max;
    }
    else
    {
        Ql[0][0] = Ql[0][0] / sqrtpi;
        Ql[1][0] = Ql[1][0] * 2.0 / sqrtpi;
        while (mp <= m_max && std::abs(Ql[mp][0]) < huge)
        {
            const auto d1 = mp + 0.5;
            Ql[mp + 1][0] = -2.0 * mp * qargu * Ql[mp][0] / d1 - (mp - 0.5) * Ql[mp - 1][0] / d1;
            ++mp;
        }
        if (mp - 1 < m_max)
            m_max = mp - 1;
        m_new = m_max;
    }
    const double fl = m_max / 2.0;
    const double cc = std::abs(double(int(fl)) - fl);
    double ar;
    if (cc < 0.4)
        ar = 1.0;
    else
        ar = -1.0;

    double gamma, gammapr;
    if (mode == 0)
    {
        gammapr = ar * gammah(m_max, huge);
        gamma = sqrtpi * gammapr;
        if (std::abs(gamma) < tiny)
        {
            std::cout << "m_max is too large for mode=0";
            std::cout << "Better try mode=1";
            return;
        }
    }
    else
    {
        gammapr = ar;
        gamma = ar;
    }

    auto dfacqs = -(gamma / Ql[m_max + 1][0]) * gamma / pi;
    // Evaluation of Pl[m_max][0],Pl[m_max+1][0]
    // we choose expansion or CF for Pl[m_max][0]
    // Depending on the values of z, m_max and mode
    int ical = 1;
    if (z / m_max > pr[iprec])
        ical = 2;
    if (z < 5.0)
        ical = 1;

    else if (z > 20.0)
    {
        if (mode != 2 && ical == 1)
            ical = 0;
    }

    if (ical == 1)
    {
        // We calculate the CF for P's
        auto fcp = fracps(qargu, m_max + 1, 0, s_eps, sqrttiny);
        const auto dd = m_max + 0.5;
        if (mode != 0)
        {
            fcp /= dd;
            dfacqs /= dd;
        }
        Pl[m_max][0] = dfacqs / qq / (1.0 - fcp * Ql[m_max][0] / Ql[m_max + 1][0]);
        Pl[m_max + 1][0] = fcp * Pl[m_max][0];
    }
    else
    {
        double Pl0;
        if (ical == 2)
            Pl0 = expan(z, mode, iprec, huge, qargu, m_max);
        else // ical == 0
            Pl0 = toroidal_asymp(z, m_max, mode, gammapr);
        Pl[m_max][0] = Pl0;
        const auto dd = m_max + 0.5;
        if (mode != 0)
            dfacqs /= dd;
        Pl[m_max + 1][0] = (Ql[m_max + 1][0] / Ql[m_max][0]) * (Pl[m_max][0] - dfacqs / qq);
    }

    // Evaluation of Pl[m_max][1]
    auto qm0 = Ql[m_max][0];
    auto dfac3 = (gamma / qm0) * gamma / pi / (0.5 - m_max);
    const auto fc = frac(z, m_max, 0, s_eps, sqrttiny);
    Pl[m_max][1] = Pl[m_max][0] * fc + dfac3;
    // Evaluation of Pl[m_max+1][1]
    qm0 = Ql[m_max + 1][0];
    double dfacc;
    if (mode == 0)
        dfacc = (0.5 + m_max);
    else
        dfacc = 1.0 / (0.5 + m_max);
    dfac3 = -(gamma / qm0) * gamma * dfacc / pi;
    const auto fc2 = frac(z, m_max + 1, 0, s_eps, sqrttiny);
    Pl[m_max + 1][1] = Pl[m_max + 1][0] * fc2 + dfac3;
    // We apply backward recurrence over m to get the set Pl[m][0],Pl[m][1]
    if (mode == 0)
    {
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][0] = -(Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / ((0.5 - mp) * (0.5 - mp));
        }
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][1] = (Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) / ((1.5 - mp) * (0.5 + mp));
        }
    }
    else
    {
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][0] = ((mp + 0.5) * Pl[mp + 1][0] + 2.0 * mp * qargu * Pl[mp][0]) / (0.5 - mp);
        }
        for (int i = 1; i <= m_max; ++i)
        {
            int mp = m_max + 1 - i;
            int m = mp - 1;
            Pl[m][1] = ((mp + 0.5) * Pl[mp + 1][1] + 2.0 * mp * qargu * Pl[mp][1]) * (mp - 0.5) / ((1.5 - mp) * (0.5 + mp));
        }
    }

    // Now, we perform the evaluation over n for each m
    // We start evaluating the set of Pl[m][n]
    for (int l = 0; l <= m_max; ++l)
    {
        int m = l;
        const double fl = m / 2.0;
        const double cc = std::abs(double(int(fl)) - fl);
        double ar;
        if (cc < 0.4)
            ar = 1.0;
        else
            ar = -1.0;

        if (mode == 0)
            gamma = gammah(m, huge) * ar * sqrtpi;
        else
            gamma = ar;

        int np = 1;
        while (np <= n_max && std::abs((np - m + 0.5) * Pl[l][np]) < huge)
        {
            Pl[l][np + 1] = (2.0 * np * z * Pl[l][np] - (np + m - 0.5) * Pl[l][np - 1]) / (np - m + 0.5);
            ++np;
        }
        n_max = np - 1;
        n_new = n_max;
    }

    // We evaluate the set of Ql[m][n]
    for (int l = 0; l <= m_max; ++l)
    {
        int m = l;
        const double fl = m / 2.0;
        const double cc = std::abs(double(int(fl)) - fl);
        if (cc < 0.4)
            ar = 1.0;
        else
            ar = -1.0;

        if (mode == 0)
            gamma = gammah(m, huge) * ar * sqrtpi;
        else
            gamma = ar;

        // We evaluate the c.f. for Q's using Lentz-Thompson
        // If m=0,1. In other case we calculate Ql[m][n_max] and
        // Ql[m][n_max+1] from the recurrence.
        if (m == 2)
        {
            if (mode == 0)
            {
                int mp = 1;
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
                int mp = 1;
                while (mp <= m_max)
                {
                    const auto d1 = mp + 0.5;
                    const auto d2 = mp - 0.5;
                    Ql[mp + 1][n_max] = -2.0 * mp * qargu * Ql[mp][n_max] / d1
                                        - (mp - n_max - 0.5) * (mp + n_max - 0.5) * Ql[mp - 1][n_max] / (d1 * d2);
                    ++mp;
                }
                mp = 1;
                while (mp <= m_max)
                {
                    const auto d1 = mp + 0.5;
                    const auto d2 = mp - 0.5;
                    Ql[mp + 1][n_max + 1] = -2.0 * mp * qargu * Ql[mp][n_max + 1] / d1
                                            - (mp - n_max - 1.5) * (mp + n_max + 0.5) * Ql[mp - 1][n_max + 1] / (d1 * d2);
                    ++mp;
                }
            }
        }
        if (m == 0 || m == 1)
        {
            const auto fc = frac(z, m, n_max, s_eps, sqrttiny);
            const auto dfac4 = factco(n_max, Pl[l][n_max + 1], m) * gamma * gamma / pi / (n_max + m + 0.5);
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
            int np = n_max + 1 - i;
            int n = np - 1;
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

    const int max_iter = 10000;
    int mm = 0;
    while (mm < max_iter)
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
    if (mm == max_iter)
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

/**
 * Evaluate @f$ P_{-1/2}^M(x) @f$ by series expansion.
 */
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

    const int max_iter = 1000;
    int k = 0;
    while (k < max_iter)
    {
        delta = (df1 + da1) * da2 * a0;
        sum += delta;
        double dcc = 0.5 + m + 2.0 * k;
        double dccp = dcc + 1.0;
        double dkk = k + 1.0;
        da2 = da2 * dccp * dcc / (dkk * dkk);
        da1 += 1.0 / dkk - 1.0 / dccp - 1.0 / dcc;
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

/**
 * Evaluation of @f$ P^m_{-1/2}(z) @f$ by asymptotic expansion.
 * @f[
 *    P^m_{-1/2}(z) \simeq (-1)^m \frac{\Gamma(m+1/2)}{\pi^{3/2}} \xi
 *     \sum_{l=0}^{1}K_l(m\alpha/2))\sum_{k=0}^{\infty}\frac{a_{2k+l}}{m^{2k+l}}
 * @f]
 * where @f$ \xi = (x-1)/2 @f$ and @f$ \alpha = ln[(x+1)/(x-1)] @f$.
 * @f$ K_0 @f$ and @f$ K_1 @f$ are the modified Bessel functions.
 */
double
toroidal_asymp(double z, int m, int mode, double gammapr)
{
    constexpr double pi = 3.14159265358979323;
    const double sqrtpi = std::sqrt(pi);
    const double psi = 0.5 * (z - 1.0);
    const double eta = (z - 1.0) / (z + 1.0);
    const double alfa = -std::log(eta);
    const double rmalfa = m * alfa;

    double sa, sb;
    toroidal_asymp_series(alfa, m, sa, sb);

    // Bessel functions K_0, K_1
    const double argu = rmalfa * 0.5;
    const auto rk0 = std::cyl_bessel_k(0.0, argu);
    const auto rk1 = std::cyl_bessel_k(1.0, argu);
    auto funci = sa * rk0 + sb * rk1;

    double gamma;
    if (mode == 0)
        gamma = gammapr / pi;
    else
        gamma = gammapr / pi / sqrtpi;
    double restof = 1.0 / std::sqrt(psi);
    double pm12 = gamma * restof * funci;
    return pm12;
}

/**
 * Inner sum for evaluation of @f$ P^m_{-1/2}(z) @f$ by asymptotic expansion.
 * @f[
 *     Sum = \sum_{k=0}^{\infty}\frac{a_{2k+1}}{m^{2k+1}}
 * @f]
 * where @f$ a_0 = \sqrt{\alpha/(e^\alpha - 1)} @f$.
 * @see GST p368 for other terms.
 */
void
toroidal_asymp_series(double a, int m, double &sa, double &sb)
{
    using Tp = double;

    // Normalization of A_0
    const Tp a0 = std::sqrt(a / (std::exp(a) - Tp{ 1 }));
    const auto a2 = a * a;

    double ak[7];

    ak[0] = Tp{ 1 };

    // ak[1] = -1/48 a + 1/2880 a^3 - 1/120960 a^5;
    ak[1] = a * (-Tp{ 1 } / Tp{ 48 } + a2 * (Tp{ 1 } / Tp{ 2880 } + a2 * (-Tp{ 1 } / Tp{ 120960 })));

    // ak[2] = 7/7680 a^2 - 13/322560 a^4;
    ak[2] = a2 * (Tp{ 7 } / Tp{ 7680 } + a2 * (-Tp{ 13 } / Tp{ 322560 }));

    // ak[3] = 7/1920 a - 571/2580480 a^3 + 1697/154828800 a^5;
    ak[3] = a * (Tp{ 7 } / Tp{ 1920 } + a2 * (-Tp{ 571 } / Tp{ 2580480 } + a2 * (Tp{ 1697 } / Tp{ 154828800 })));

    // ak[4] = -31/64512 a^2 + 22763/495452160 a^4;
    ak[4] = a2 * (-Tp{ 31 } / Tp{ 64512 } + a2 * (Tp{ 22763 } / Tp{ 495452160 }));

    // ak[5] = -31/16128 a + Tp{7691}/Tp{30965760} a^3 - Tp{5501381} / Tp{261598740480} a^5;
    ak[5] = a * (-Tp{ 31 } / Tp{ 16128 } + a2 * (Tp{ 7691 } / Tp{ 30965760 } + a2 * (-Tp{ 5501381 } / Tp{ 261598740480 })));

    // ak[6] = 127/245760 a^2 - 44593/519045120 a^4;
    ak[6] = a2 * (Tp{ 127 } / Tp{ 245760 } + a2 * (-Tp{ 44593 } / Tp{ 519045120 }));

    auto rmk = Tp{ 1 };
    sa = Tp{ 0 };
    sb = Tp{ 0 };
    for (int k = 0; k <= 6; ++k)
    {
        const auto fl = Tp(k) / Tp{ 2 };
        const auto cc = std::abs(double(int(fl)) - fl);
        if (cc < 0.4)
            sa += ak[k] / rmk;
        else
            sb += ak[k] / rmk;
        rmk *= m;
    }
    sa *= a0;
    sb *= a0;
    return;
}

int
main()
{
    unsigned int m_max = 10;
    unsigned int n_max = 10;
    std::vector<std::vector<double>> P(m_max + 2, std::vector<double>(n_max + 2));
    std::vector<std::vector<double>> Q(m_max + 2, std::vector<double>(n_max + 2));
    int m_new, n_new;
    for (double x : {1.01, 2.0, 5.0})
    {
        toroidal_harmonic(x, m_max, n_max, P, Q, m_new, n_new);
        std::cout << "\nx = " << x;
        for (unsigned int m = 0; m <= m_max; ++m)
        {
            std::cout << '\n';
            for (unsigned int n = 0; n <= n_max; ++n)
            {
                std::cout << ' ' << std::setw(12) << P[m][n];
            }
        }
    }
    std::cout << '\n';
}

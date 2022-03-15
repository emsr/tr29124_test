#ifndef BRENT_TCC
#define BRENT_TCC 1

#include <cmath>
#include <numbers>
#include <algorithm> // minmax
#include <stdexcept>

namespace emsr
{

/**
 * 
 */
template<typename Func, typename Real>
FuncPt<Real>
brent(Func func, Real ax, Real bx, Real cx, Real tol)
{
    constexpr auto gold = std::numbers::phi_v<Real>;
    constexpr auto gr = gold - Real{1};
    constexpr auto gc = Real{1} - gr;
    constexpr auto zeps = Real{1.0e-10};
    constexpr auto zero = Real{0};
    constexpr auto half = Real{0.5};
    const int max_iter = 100;

    auto a = std::min(ax, cx);
    auto b = std::max(ax, cx);
    auto x = FuncPt(func, bx);
    auto w = x;
    auto v = x;
    auto e = zero;
    for (int iter = 1; iter <= max_iter; ++iter)
    {
        const auto x_mid = half * (a + b);
        const auto tol1 = tol * std::abs(x.arg) + zeps;
        const auto tol2 = Real{2} * tol1;

        if (std::abs(x.arg - x_mid) <= tol2 - half * (b - a))
        {
            return x;
        }

        Real d;
        if (std::abs(e) > tol1)
        {
            auto r = (x.arg - w.arg) * (x.val - v.val);
            auto q = (x.arg - v.arg) * (x.val - w.val);
            auto p = (x.arg - v.arg) * q - (x.arg - w.arg) * r;
            q = Real{2} * (q - r);
            if (q > zero)
                p = -p;
            q = std::abs(q);
            auto e_temp = e;
            e = d;
            if (std::abs(p) >= std::abs(half * q * e_temp)
                || p <= q * (a - x.arg) || p >= q * (b - x.arg))
            {
                e = x.arg >= x_mid ? a - x.arg : b - x.arg;
                d = gc * e;
            }
            else
            {
                d = p / q;
                auto uarg = x.arg + d;
                if (uarg - a < tol2 || b - uarg < tol2)
                    d = std::copysign(tol1, x_mid - x.arg);
            }
        }
        else
        {
            e = x.arg >= x_mid ? a - x.arg : b - x.arg;
            d = gc * e;
        }

        const auto uarg = std::abs(d) >= tol1
                        ? x.arg + d
                        : x.arg + std::copysign(tol1, d);
        auto u = FuncPt(func, uarg);
        if (u.val < x.val)
        {
            if (u.arg >= x.arg)
                a = x.arg;
            else
                b = x.arg;
            v = w;
            w = x;
            x = u;
        }
        else
        {
            if (u.arg < x.arg)
                a = u.arg;
            else
                b = u.arg;
            if (u.val <= w.val || w.arg == x.arg)
            {
                v = w;
                w = u;
            }
            else if (u.val <= v.val || v.arg == x.arg || v.arg == w.arg)
            {
                v = u;
            }
        }
    }
    throw std::logic_error("brent: Maximum iterations exceeded.");
}

/**
 * 
 */
template<typename Func, typename Deriv, typename Real>
DerivPt<Real>
brent(Func func, Deriv deriv, Real ax, Real bx, Real cx, Real tol)
{
    constexpr auto zeps = Real{1.0e-10};
    constexpr auto zero = Real{0};
    constexpr auto half = Real{0.5};
    const int max_iter = 100;

    auto a = std::min(ax, cx);
    auto b = std::max(ax, cx);
    auto x = DerivPt(func, deriv, bx);
    auto w = x;
    auto v = x;
    auto e = zero;
    for (int iter = 1; iter <= max_iter; ++iter)
    {
        const auto x_mid = half * (a + b);
        const auto tol1 = tol * std::abs(x.arg) + zeps;
        const auto tol2 = Real{2} * tol1;

        if (std::abs(x.arg - x_mid) <= tol2 - half * (b - a))
        {
            return x;
        }

        Real d;
        if (std::abs(e) > tol1)
        {
            auto d1 = Real{2} * (b - a);
            auto d2 = d1;
            if (w.der != x.der)
                d1 = (w.arg - x.arg) * x.der / (x.der - w.der);
            if (v.der != x.der)
                d2 = (v.arg - x.arg) * x.der / (x.der - v.der);
            auto u1 = x.arg + d1;
            auto u2 = x.arg + d2;
            auto ok1 = (a - u1) * (u1 - b) > zero && x.der * d1 <= zero;
            auto ok2 = (a - u2) * (u2 - b) > zero && x.der * d2 <= zero;
            const auto e_temp = e;
            e = d;
            if (ok1 || ok2)
            {
                if (ok1 && ok2)
                    d = std::abs(d1) <= std::abs(d2) ? d1 : d2;
                else if (ok1)
                    d = d1;
                else
                    d = d2;

                if (std::abs(d) <= std::abs(half * e_temp))
                {
                    auto uarg = x.arg + d;
                    if (uarg - a < tol2 || b - uarg < tol2)
                        d = std::copysign(tol1, x_mid - x.arg);
                }
                else
                {
                    e = x.der >= zero ? a - x.arg : b - x.arg;
                    d = half * e;
                }
            }
            else
            {
                e = x.der >= zero ? a - x.arg : b - x.arg;
                d = half * e;
            }
        }
        else
        {
            e = x.der >= zero ? a - x.arg : b - x.arg;
            d = half * e;
        }

        const auto uarg = std::abs(d) >= tol1
                        ? x.arg + d
                        : x.arg + std::copysign(tol1, d);
        DerivPt<Real> u(func, deriv, uarg);
        if (std::abs(d) < tol1 && u.val > x.val)
        {
            return x;
        }

        if (u.val <= x.val)
        {
            if (u.arg >= x.arg)
                a = x.arg;
            else
                b = x.arg;
            v = w;
            w = x;
            x = u;
        }
        else
        {
            if (u.arg < x.arg)
                a = u.arg;
            else
                b = u.arg;
            if (u.val <= w.val || w.arg == x.arg)
            {
                v = w;
                w = u;
            }
            else if (u.val <= v.val || v.arg == x.arg || v.arg == w.arg)
            {
                v = u;
            }
        }
    }
    throw std::logic_error("brent: Maximum iterations exceeded.");
}

} // namespace emsr

#endif // BRENT_TCC

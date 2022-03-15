#ifndef GOLDEN_TCC
#define GOLDEN_TCC 1

#include <cmath>
#include <numbers>

namespace emsr
{

/**
 * 
 */
template<typename Func, typename Real>
void
mean_bracket(Func func, MinBracket<Real>& mb)
{
    constexpr auto gold = std::numbers::phi_v<Real>;
    constexpr auto tiny = 100 * std::numeric_limits<Real>::min();
    const auto glimit = Real{100};

    // We want to go downhill from a to b.
    mb.a = FuncPt(func, mb.a.arg);
    mb.b = FuncPt(func, mb.b.arg);
    if (mb.b.val > mb.a.val)
        std::swap(mb.a, mb.b);
    mb.c = FuncPt(func, mb.b.arg + gold * (mb.b.arg - mb.a.arg));
    while (mb.b.val > mb.c.val)
    {
        // Parabolic interpolation.
        const auto r = (mb.b.arg - mb.a.arg) * (mb.b.val - mb.c.val);
        const auto q = (mb.b.arg - mb.c.arg) * (mb.b.val - mb.a.val);
        const auto qmr_safe = std::copysign(std::max(std::abs(q - r), tiny), q - r);
        const auto ulim = mb.b.arg * glimit * (mb.c.arg - mb.b.arg);
        const auto uarg = mb.b.arg
                        - ((mb.b.arg - mb.c.arg) * q - (mb.b.arg - mb.a.arg) * r)
                          / (Real{2} * qmr_safe);
        auto u = FuncPt(func, uarg);
        if ((mb.b.arg - u.arg) * (u.arg - mb.c.arg) > Real{0})
        {
            if (u.val < mb.c.val)
            {
                mb.a = mb.b;
                mb.b = u;
                return;
            }
            else if (u.val > mb.b.val)
            {
                mb.c = u;
                return;
            }
        }
        else if ((mb.c.arg - u.arg) * (u.arg - ulim) > Real{0})
        {
            if (u.val < mb.c.val)
            {
                mb.b = mb.c;
                mb.c = u;
                u = FuncPt(func, mb.c.arg + gold * (mb.c.arg - mb.b.arg));
            }
        }
        else if ((u.arg - ulim) * (ulim - mb.c.arg) > Real{0})
        {
            u = FuncPt(func, ulim);
        }
        else
        {
            u = FuncPt(func, mb.c.arg + gold * (mb.c.arg - mb.b.arg));
        }

        mb.a = mb.b;
        mb.b = mb.c;
        mb.c = u;
    }
}

/**
 * 
 */
template<typename Func, typename Real>
FuncPt<Real>
golden(Func func, MinBracket<Real>& mb, Real tol)
{
    constexpr auto gold = std::numbers::phi_v<Real>;
    constexpr auto gr = gold - Real{1};
    constexpr auto gc = Real{1} - gr;
    constexpr int max_iter = 100;

    auto x_min = mb.a.arg;
    auto x_max = mb.c.arg;
    Real x1, x2;
    if (std::abs(mb.c.arg - mb.b.arg) > std::abs(mb.b.arg - mb.a.arg))
    {
        x1 = mb.b.arg;
        x2 = mb.b.arg + gc * (mb.c.arg - mb.b.arg);
    }
    else
    {
        x2 = mb.b.arg;
        x1 = mb.b.arg - gc * (mb.b.arg - mb.a.arg);
    }
    auto pt1 = FuncPt(func, x1);
    auto pt2 = FuncPt(func, x2);
    int iter = 0;
    while (std::abs(x_max - x_min) > tol * (std::abs(pt1.arg) + std::abs(pt2.arg)))
    {
        ++iter;
        if (iter > max_iter)
            throw std::logic_error("golden: Maximum iterations exceeded.");

        if (pt2.val < pt1.val)
        {
            x_min = pt1.arg;
            pt1 = pt2;
            pt2 = FuncPt(func, gr * pt1.arg + gc * x_max);
        }
        else
        {
            x_max = pt2.arg;
            pt2 = pt1;
            pt1 = FuncPt(func, gr * pt2.arg + gc * x_min);
        }
    }
    if (pt1.val < pt2.val)
    {
        return pt1;
    }
    else
    {
        return pt2;
    }
}

} // namespace emsr

#endif // GOLDEN_TCC

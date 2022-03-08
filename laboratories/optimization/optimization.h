#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H 1

namespace emsr
{

/**
 * A struct for a function argument and value.
 */
template<typename Real>
struct FuncPt
{
    Real arg;
    Real val;

    FuncPt() = default;

    FuncPt(Real x, Real f)
    : arg(x), val(f)
    { }

    template<typename Func>
    FuncPt(Func func, Real x)
    : arg(x), val(func(x))
    { }
};

/**
 * A struct for the function minimum mb.
 */
template<typename Real>
struct MinBracket
{
    FuncPt<Real> a;
    FuncPt<Real> b;
    FuncPt<Real> c;

    MinBracket() = default;

    // Start the first two points.
    template<typename Func>
    MinBracket(Func func, Real a0, Real b0)
    : a(func, a0), b(func, b0), c()
    { }

    // Try all three points.
    template<typename Func>
    MinBracket(Func func, Real a0, Real b0, Real c0)
    : a(func, a0), b(func, b0), c(func, c0)
    { }
};

/**
 * A struct for a function argument, value and derivative.
 */
template<typename Real>
struct DerivPt
{
    Real arg;
    Real val;
    Real der;

    DerivPt() = default;

    DerivPt(Real x, Real f, Real d)
    : arg(x), val(f), der(d)
    { }

    template<typename Func, typename Deriv>
    DerivPt(Func func, Deriv deriv, Real x)
    : arg(x), val(func(x)), der(deriv(x))
    { }
};

}

#include <golden.h>

#include <brent.h>

#endif // OPTIMIZATION_H

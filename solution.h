#ifndef SOLUTION_H
#define SOLUTION_H 1

#include <complex>
#include <variant>
#include <iosfwd>

namespace __gnu_cxx
{

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

  template<typename _Real>
    constexpr bool
    is_valid(const solution_t<_Real>& __x)
    { return __x.index() != 0; }

  template<typename _Real>
    constexpr _Real
    real(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return std::numeric_limits<_Real>::quiet_NaN();
      else if (__x.index() == 1)
	return std::get<1>(__x);
      else
	return std::real(std::get<2>(__x));
    }

  template<typename _Real>
    constexpr _Real
    imag(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return std::numeric_limits<_Real>::quiet_NaN();
      else if (__x.index() == 1)
	return _Real{0};
      else
	return std::imag(std::get<2>(__x));
    }

  /**
   * Return the solution as a complex number or NaN.
   */
  template<typename _Real>
    constexpr std::complex<_Real>
    to_complex(const solution_t<_Real>& __x)
    {
      constexpr auto _S_NaN = std::numeric_limits<_Real>::quiet_NaN();
      if (__x.index() == 0)
	return std::complex<_Real>{_S_NaN, _S_NaN};
      else if (__x.index() == 1)
	return std::complex<_Real>(std::get<1>(__x), _Real{0});
      else
	return std::get<2>(__x);
    }

  /**
   * Addition operators...
   */
  template<typename _Real>
    constexpr std::complex<_Real>
    operator+(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      constexpr auto _S_NaN = std::numeric_limits<_Real>::quiet_NaN();
      if (__x.index() == 0 || __y.index() == 0)
	return std::complex<_Real>{_S_NaN, _S_NaN};
      else
	return to_complex(__x) + to_complex(__y);
    }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator+(const solution_t<_Real>& __x, _Real __y)
    { return operator+(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator+(_Real __x, const solution_t<_Real>& __y)
    { return operator+(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator+(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator+(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator+(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator+(solution_t<_Real>(__x), __y); }

  /**
   * Subtraction operators...
   */
  template<typename _Real>
    constexpr std::complex<_Real>
    operator-(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      constexpr auto _S_NaN = std::numeric_limits<_Real>::quiet_NaN();
      if (__x.index() == 0 || __y.index() == 0)
	return std::complex<_Real>{_S_NaN, _S_NaN};
      else
	return to_complex(__x) - to_complex(__y);
    }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator-(const solution_t<_Real>& __x, _Real __y)
    { return operator-(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator-(_Real __x, const solution_t<_Real>& __y)
    { return operator-(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator-(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator-(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator-(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator-(solution_t<_Real>(__x), __y); }

  /**
   * Multiplication operators...
   */
  template<typename _Real>
    constexpr std::complex<_Real>
    operator*(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      constexpr auto _S_NaN = std::numeric_limits<_Real>::quiet_NaN();
      if (__x.index() == 0 || __y.index() == 0)
	return std::complex<_Real>{_S_NaN, _S_NaN};
      else
	return to_complex(__x) * to_complex(__y);
    }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator*(const solution_t<_Real>& __x, _Real __y)
    { return operator*(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator*(_Real __x, const solution_t<_Real>& __y)
    { return operator*(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator*(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator*(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator*(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator*(solution_t<_Real>(__x), __y); }

  /**
   * division operators...
   */
  template<typename _Real>
    constexpr std::complex<_Real>
    operator/(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      constexpr auto _S_NaN = std::numeric_limits<_Real>::quiet_NaN();
      if (__x.index() == 0 || __y.index() == 0)
	return std::complex<_Real>{_S_NaN, _S_NaN};
      else
	return to_complex(__x) / to_complex(__y);
    }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator/(const solution_t<_Real>& __x, _Real __y)
    { return operator/(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator/(_Real __x, const solution_t<_Real>& __y)
    { return operator/(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator/(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator/(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr std::complex<_Real>
    operator/(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator/(solution_t<_Real>(__x), __y); }

  /**
   * Test for equality and inequality.
   */
  template<typename _Real>
    constexpr bool
    operator==(const solution_t<_Real>& __x,
	       const solution_t<_Real>& __y)
    {
      if (__x.index() == 0 || __y.index() == 0)
	return false;
      else if (__x.index() == __y.index())
	{
	  if (__x.index() == 1)
	    return std::get<1>(__x) == std::get<1>(__y);
	  else
	    return std::get<2>(__x) == std::get<2>(__y);
	}
      else
	return false;
    }

  template<typename _Real>
    bool
    operator==(const solution_t<_Real>& __x, _Real __y)
    { return __x == solution_t<_Real>(__y); }

  template<typename _Real>
    constexpr bool
    operator==(_Real __x, const solution_t<_Real>& __y)
    { return solution_t<_Real>(__x) == __y; }

  template<typename _Real>
    constexpr bool
    operator==(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return __x == solution_t<_Real>(__y); }

  template<typename _Real>
    constexpr bool
    operator==(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return solution_t<_Real>(__x) == __y; }

  template<typename _Real>
    constexpr bool
    operator!=(const solution_t<_Real>& __x,
	       const solution_t<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    bool
    operator!=(const solution_t<_Real>& __x, _Real __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(_Real __x, const solution_t<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return !(__x == __y); }

  /**
   * Lexicographic order of solutions as complex numbers.
   * Null solutions compare as less than except to another null solution.
   *
   * A tribool might be a good thing for this when either
   * of the solutions is null.
   */
  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0 && __y.index() == 0)
	return false;
      else if (__x.index() == 0)
	return true;
      else if (__y.index() == 0)
	return false;
      else
	{
	  const auto __rex = real(__x);
	  const auto __rey = real(__y);
	  if (std::get<1>(__rex) < std::get<1>(__rey))
	    return true;
	  else if (__rex == __rey)
	    return std::get<1>(imag(__x)) < std::get<1>(imag(__y));
	  else
	    return false;
	}
    }

  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, _Real __y)
    { return operator<(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr bool
    operator<(_Real __x, const solution_t<_Real>& __y)
    { return operator<(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return operator<(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr bool
    operator<(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator<(solution_t<_Real>(__x), __y); }

} // namespace __gnu_cxx

  /**
   * Output a solution to a stream.
   */
  template<typename _Char, typename _Real>
    std::basic_ostream<_Char>&
    operator<<(std::basic_ostream<_Char>& __out,
	       const __gnu_cxx::solution_t<_Real>& __sln)
    {
      const auto __idx = __sln.index();
      if (__idx == 0)
	__out << "null";
      else if (__idx == 1)
	__out << std::get<1>(__sln);
      else if (__idx == 2)
	__out << std::get<2>(__sln);
      return __out;
    }

#endif // SOLUTION_H

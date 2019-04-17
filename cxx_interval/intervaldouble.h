#ifndef INC_INTERVAL
#define INC_INTERVAL

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2015
 *                       Future Team Aps 
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Future Team Software License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :   intervaldouble.h
 * Module ID Nbr   :   
 * Description     :   Interval arithmetic template class
 *                     Works with both float and double
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/28dec14		Initial release	
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/


/* define version string */
static char _VinterP_[] = "@(#)intervaldouble.h 01.01 -- Copyright (C) Future Team Aps";

#include <float.h>
#include <algorithm>

// By default HARDWARE_SUPPORT controlled if IEEE754 floating point control can be used for interval arithmetic.
// The intervaldouble.h requires this to be defined.
#define HARDWARE_SUPPORT

/// The four different interval classification
/// # ZERO			a=0 && b=0
/// # POSITIVE		a>=0 && b>0
/// # NEGATIVE		a<0 && b<=0
/// # MIXED			a<0 && b>0
enum int_class { NO_CLASS, ZERO, POSITIVE, NEGATIVE, MIXED };

/// The four different ronding modes
/// # ROUND_NEAR  Rounded result is the closest to the infinitely precise result.
/// # ROUND_DOWN  Rounded result is close to but no greater than the infinitely precise result.
/// # ROUND_UP    Rounded result is close to but no less than the infinitely precise result.
/// # ROUND ZERO  Rounded result is close to but no greater in absolute value than the infinitely precise result.
enum round_mode { ROUND_NEAR, ROUND_UP, ROUND_DOWN, ROUND_ZERO };

// 
// Interval class
// Realistically the class Type can be float, double. Any other type is not supported
// Since float and double are done unsing the Intel cpu (H/W) and using "specilization" . 
//
template<class _IT> class interval {
   _IT low, high;
   public:
      typedef _IT value_type;

      // constructor. zero, one or two arguments for type _IT
      interval()							{ low = _IT(0); high = _IT(0); }
	  interval( const _IT& d )				{ low = _IT(d); high = _IT(d); }
	  interval( const _IT& l, const _IT& h) { if( l < h ) { low =l; high = h; } else { low = h; high = l; } }
	  // Constrcutor for mixed type _IT != _X (base types). Allows auto construction of e.g. interval<float_precision> x(float)
	 template <class _X> interval(const _X& x) { low = _IT(x); high = _IT(x); }

      // constructor for any other type to _IT. Both up and down conversion possible
      template<class X> interval( const interval<X>& a ) 
		{
		if (a.lower() < a.upper()) { fpdown();  low = _IT(a.lower()); fpup();  high = _IT(a.upper()); fpnear();  }
		else { fpdown();  low = _IT(a.upper()); fpup();  high = _IT(a.lower()); fpnear(); }
		}
	 
      // Coordinate functions
      _IT upper() const					{ return high; }
      _IT lower() const					{ return low; }
      _IT upper( const _IT& u )			{ return ( high = u ); }
      _IT lower( const _IT& l )			{ return ( low = l ); }

      _IT center() const				{ return ( high + low ) / _IT(2); }
      _IT radius() const				{ _IT r; r =( high - low ) / _IT(2); if( r < _IT(0) ) r = -r; return r; }
	  _IT width() const					{ _IT r; r = high - low; if (r < _IT(0)) r = -r; return r; }
	
	  bool contain_zero() const			{ return low <= _IT(0) && _IT(0) <= high; }  // Obsolete. use contains() instead.
	  bool contain( const _IT& f=_IT(0)){ return low <= f && f <= high;  }
	  bool contain(const interval<_IT>& i) { return low <= i.lower() && i.upper() <= high; }
	  bool is_empty() const				{ return high < low; }
	
	  enum int_class is_class() const	{ 
										if (low == _IT(0) && high == _IT(0)) return ZERO;
										if (low >= _IT(0) && high > _IT(0)) return POSITIVE;
										if (low < _IT(0) && high <= _IT(0)) return NEGATIVE;
										if (low < _IT(0) && high > _IT(0)) return MIXED;
										return NO_CLASS;
										}

	  // Operators
	  operator short() const			{ return (short)((high + low) / _IT(2)); }				// Conversion to short 	 
	  operator int() const				{ return (int)center() /*(( high + low ) / _IT(2) )*/; }			// Conversion to int 
	  operator long() const				{ return (long)((high + low) / _IT(2)); }				// Conversion to long 
	  operator unsigned short() const	{ return (unsigned short)((high + low) / _IT(2)); }		// Conversion to unsigned short 	 
	  operator unsigned int() const		{ return (unsigned int)((high + low) / _IT(2)); }		// Conversion to unsigned int 
	  operator unsigned long() const	{ return (unsigned long)((high + low) / _IT(2)); }		// Conversion to unsigned long 
	  operator double() const			{ return (double)( ( high + low ) / _IT(2) ); }			// Conversion to double 
	  operator float() const			{ return high == low? (float)low : (float)((high + low) / _IT(2)); }				// Conversion to float 

      _IT *ref_lower()					{ return &low; }
      _IT *ref_upper()					{ return &high; }

      // Essential operators
      interval<_IT>& operator= ( const interval<_IT>& );
      interval<_IT>& operator+=( const interval<_IT>& );
      interval<_IT>& operator-=( const interval<_IT>& );
      interval<_IT>& operator*=( const interval<_IT>& );
      interval<_IT>& operator/=( const interval<_IT>& );
	  interval<_IT>& operator&=( const interval<_IT>& );
	  interval<_IT>& operator|=( const interval<_IT>& );
	  interval<_IT>& operator^=( const interval<_IT>& );
	  
	  // Exception class. No used
	  class bad_int_syntax {};
	  class bad_float_syntax {};
	  class out_of_range   {};
	  class divide_by_zero {};
	  class domain_error   {};
	  class base_error		{};
   };


// Unary and Binary arithmetic
// Arithmetic + Binary and Unary
template <class _IT, class _X> interval<_IT> operator+( const interval<_IT>&, const _X&);
template <class _IT, class _X> interval<_IT> operator+( const _X&, const interval<_IT>&);
template<class _IT> interval<_IT> operator+( const interval<_IT>&, const interval<_IT>& ); 
template<class _IT> interval<_IT> operator+( const interval<_IT>& );									// Unary 

// Arithmetic - Binary and Unary
template <class _IT, class _X> interval<_IT> operator-(const interval<_IT>&, const _X&);
template <class _IT, class _X> interval<_IT> operator-(const _X&, const interval<_IT>&);
template<class _IT> interval<_IT> operator-( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator-( const interval<_IT>& );									// Unary

// Arithmetic * Binary
template <class _IT, class _X> interval<_IT> operator*(const interval<_IT>&, const _X&);
template <class _IT, class _X> interval<_IT> operator*(const _X&, const interval<_IT>&);
template<class _IT> interval<_IT> operator*( const interval<_IT>&, const interval<_IT>& );

// Arithmetic / Binary
template <class _IT, class _X> interval<_IT> operator/(const interval<_IT>&, const _X&);
template <class _IT, class _X> interval<_IT> operator/(const _X&, const interval<_IT>&);
template<class _IT> interval<_IT> operator/( const interval<_IT>&, const interval<_IT>& );

// Boolean Comparison Operators
template <class _IT, class _X> bool operator==(const interval<_IT>&, const _X&);
template <class _IT, class _X> bool operator==(const _X&, const interval<_IT>&);
template <class _IT, class _X> bool operator==(const interval<_IT>&, const interval<_X>&);
//template<class _IT> bool operator==(const interval<_IT>&, const interval<_IT>&);

template <class _IT, class _X> bool operator!=(const interval<_IT>&, const _X&);
template <class _IT, class _X> bool operator!=(const _X&, const interval<_IT>&);
template <class _IT, class _X> bool operator!=(const interval<_IT>&, const interval<_X>&);
//template<class _IT> bool operator!=(const interval<_IT>&, const interval<_IT>&);

// Other functions
template<class _IT> interval<_IT> abs(const interval<_IT>&);

// Manifest Constants like PI, LN2 and LN10
inline interval<float> int_pifloat();
inline interval<double> int_pidouble();
inline interval<float> int_ln2float();
inline interval<double> int_ln2double();
inline interval<float> int_ln10float();
inline interval<double> int_ln10double();

// Elementary functions
inline interval<float> sqrt(const interval<float>&);
inline interval<double> sqrt(const interval<double>&);
inline interval<float> log( const interval<float>& );
inline interval<double> log(const interval<double>&);
inline interval<float> log10(const interval<float>&);
inline interval<double> log10(const interval<double>&);
inline interval<float> exp(const interval<float>&, const float);
inline interval<double> exp(const interval<double>&, const double);
inline interval<float> pow(const interval<float>&, const float );
inline interval<double> pow( const interval<double>&, const double );

// Trigonometric functions
inline interval<float> sin(const interval<float>&);
inline interval<double> sin(const interval<double>&);
inline interval<float> cos(const interval<float>&);
inline interval<double> cos(const interval<double>&);
inline interval<float> tan(const interval<float>&);
inline interval<double> tan(const interval<double>&);
inline interval<float> asin(const interval<float>&);
inline interval<double> asin(const interval<double>&);
inline interval<float> acos(const interval<float>&);
inline interval<double> acos(const interval<double>&);
inline interval<float> atan(const interval<float>&);
inline interval<double> atan(const interval<double>&);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//	Floating point control for the IEEE754 hardware. Only fo non managed application
//	Enable by defined #define HARDWARE_SUPPORT
//  Currently it only works with HARDWARE_SUPPORT defined.
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef HARDWARE_SUPPORT
inline void fpnear()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_NEAR, _MCW_RC);
	cc=__control87_2(_RC_NEAR, _MCW_RC, &f87_cw, &sse2_cw );
	}

inline void fpdown()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_DOWN, _MCW_RC);
	cc=__control87_2(_RC_DOWN, _MCW_RC, &f87_cw, &sse2_cw);
	}

inline void fpup()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_UP, _MCW_RC);
	cc=__control87_2(_RC_UP, _MCW_RC, &f87_cw, &sse2_cw);
	}
#else
inline void fpnear()	{}
inline void fpdown()	{}
inline void fpup()		{}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//   End Floating point control for the IEEE754 hardware
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Output Operator <<
//
template<class _Ty> inline std::ostream& operator<<( std::ostream& strm, interval<_Ty>& a )
	{ return strm << "[" << a.lower() << "," << a.upper() << "]"; }

// Input operator >>
//
template<class _Ty> inline std::istream& operator>>( std::istream& strm, interval<_Ty>& c ) 
   {
   _Ty l, u; char ch;
   if( strm >> ch && ch != '[')
      strm.putback(ch), strm >> l, u = l;
	else
      if( strm >> l >> ch && ch != ',')
	      if( ch == ']')
	         u = l;
	      else 
            strm.putback( ch ); // strm.setstate(std::ios::failbit);
	   else
         if( strm >> u >> ch && ch != ']')
	         strm.putback( ch ); //, strm.setstate(ios_base::failbit);
	
   if(!strm.fail())
	   c = interval<_Ty>( l, u );

   return strm;
   }
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Essential Operators =,+=,-=,*=,/=
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


// Assignment operator. Works for all class types
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator=( const interval<_IT>& a )
   {
   low = a.lower();
   high = a.upper();
   return *this;
   }

// += operator. Works all other classes. 
// Please note that this is for all integer classes. interval<int>, interval<long>
// were there os no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator+=( const interval<_IT>& a )
   {
   fpdown();
   low += a.lower();
   fpup();
   high += a.upper();
   fpnear();
   return *this;
   }

// -= operator. Works all other classes. 
// Please note that this is for all integer classes. interval<int>, interval<long>
// were there is no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator-=( const interval<_IT>& a )
   {
	fpdown();
   low -= a.high;
   fpup();
   high -= a.low;
   fpnear();
   return *this;
   }

// Works all other classes. 
// Please note that this is for all interger classes. interval<int>, interval<long>, 
// were there is no loss of precision
// Instead of doing the mindless low = MIN(low*a.high, low*a.low,high*a.low,high*a.high) and
// high = MAX(low*a.high, low*a.low,high*a.low,high*a.high) requiring a total of 8 multiplication
//
//   low, high, a.low, a.high    result
//    +     +     +     +        +  +  [ low*a.low, high*a.high ]
//    +     +     -     +        -  +  [ high*a.low, high*a.high ]
//    +     +     -     -        -  -  [ high*a.low, low*a.high ]
//    -     +     +     +        -  +  [ low*a.high, high*a.high ]  
//    -     +     -     +        -  +  [ MIN(low*a.high,high*a.low), MAX(low*a.low,high*a.high) ]
//    -     +     -     -        -  -  [ high*a.low, low*a.low ]
//    -     -     +     +        -  -  [ low*a.high, high,a.low ]
//    -     -     -     +        -  -  [ low*a.high, low*a.low ]
//    -     -     -     -        +  +  [ high*a.high, low * a.low ]
//
 template<class _IT> inline interval<_IT>& interval<_IT>::operator*=( const interval<_IT>& a )
	{
	_IT l, h, t;

	if (low >= 0) // 
	{ // both low and high >= 0
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		{ // a.low and a.high < 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.upper();
		}
	}
	else
	if (high >= 0)
	{  // low < 0, high >= 0
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper(); if (l > (t = high * a.lower())) l = t;
			fpup();
			h = high * a.upper(); if (h < (t = low * a.lower())) h = t;
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.lower();
		}
	}
	else
	{ // low and high are < 0 
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.lower();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = low * a.lower();
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.upper();
			fpup();
			h = low * a.lower();
		}
	}

	low = l;
	high = h;
	fpnear();

	return *this;
	}

// Works for all other classes
// Please note that this is for all interger classes. interval<int>, interval<long>
// were there is no loss of precision
// Actually there is specialization for both <int>
template<class _IT> inline interval<_IT>& interval<_IT>::operator/=( const interval<_IT>& b )
   {
   interval<_IT> a, c;

   fpdown();
   c.low = (_IT)1 / b.upper();
   fpup();
   c.high = (_IT)1 / b.lower();
   fpnear();
   a = interval( low, high );
   c *= a;

   low = c.lower();
   high = c.upper();

   return *this; 
   }

// Specialization for int and /=
//
inline interval<int>& interval<int>::operator/=( const interval<int>& b )
   {
   double tlow, thigh;
   interval<int> a;
   interval<double> c;

   tlow = 1 / (double)b.upper();
   thigh = 1 / (double)b.lower();

   a = interval( low, high );
   c = interval<double>( tlow, thigh );
   c *= a;

   low = (int)floor( c.lower() );
   high = (int)ceil( c.upper() );

   return *this; 
   }

// Works on all classes. 
// Return the intersection
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator&=(const interval<_IT>& a)
	{
	if (a.lower() > low ) 
		low = a.lower();
	if (a.upper() < high)
		high = a.upper();
	if (low > high)  // Empty set
		{
		low = 0; high = 0;
		}

	return *this;
	}

// Works on all classes. 
// Return the union
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator|=(const interval<_IT>& a)
	{
	if (low > a.upper() || high < a.lower())
		{
		if (a.upper() - a.lower() > high - low)
			{ // return the largest set
			low = a.lower();
			high = a.upper();
			}
		}
	else
		{ // non empty intersection
		if (a.lower() < low)
			low = a.lower();
		if (a.upper() > high)
			high = a.upper();
		}
	}

	
// Works on all classes. 
// Return the set minus
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator^=(const interval<_IT>& a)
	{
	if ( a.lower() < high && a.upper() > low ) // intersection is not empty
		{
		if (a.upper() <= low)
			low = a.upper();
		else
			if (a.high() >= high)
				high = a.lower();
		}

	return *this;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Essential Operators
///
//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Binary and Unary Operators +,-,*,/
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c += interval<_IT>(_IT(b));
	return c; 
	}

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+( const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c += interval<_IT>(_IT(a));
	return c;
	}

// Binary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c += b;
	return c;
	}


// Unary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+( const interval<_IT>& a )
   {
   return a;
   }

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c -= interval<_IT>(_IT(b));
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c -= b;
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c -= b;
   return c;
   }


// Unary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a )
   {
   interval<_IT> c(0);

   c -= a;
   return c;
   }


// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c *= interval<_IT>(_IT(b));
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c *= interval<_IT>(_IT(a));
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator*( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c *= b;
   return c;
   }

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c /= interval<_IT>(_IT(b));
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c /= b;
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator/( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   if ( c == b && b.is_class() != ZERO )
	  c = interval<_IT>(1,1);
   else
      c /= b;
   
   return c;
   }

// Binary & operator
// Return intersection
// Works for all classes
//
template<class _IT> inline interval<_IT> operator&( const interval<_IT>& a, const interval<_IT>& b )
	{
	interval<_IT> c(a);

	c &= b;
	return c;
	}

// Binary | operator.
// Return union
// Works for all classes
//
template<class _IT> inline interval<_IT> operator|(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c |= b;
	return c;
	}

// Binary ^ operator
// Return set minus
// Works for all classes
//
template<class _IT> inline interval<_IT> operator^(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c ^= b;
	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Binary and Unary Operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Boolean Interval for == and !=
///
//////////////////////////////////////////////////////////////////////////////////////


// Binary == operator
// Works for all mixed classes
//
template<class _IT, class _X> inline bool operator==(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(b);
	return c.lower() == a.lower() && c.upper() == a.upper();
	}

// Binary == operator
// Works for all mixed classes
//
template<class _IT, class _X> inline bool operator==(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);
	return c.lower() == b.lower() && c.upper() == b.upper();
	}

// == operator
// Works for all classes
//
template<class _IT, class _X> inline bool operator==(const interval<_IT>& a, const interval<_X>& b)
	{
	return a.lower() == b.lower() && a.upper() == b.upper();
	}

// == operator
// Works for all classes
//
//template<class _IT> inline bool operator==(const interval<_IT>& a, const interval<_IT>& b)
//{
//	return a.lower() == b.lower() && a.upper() == b.upper();
//}


// Binary != operator
// Works for all mixed classes
//
template<class _IT, class _X> inline bool operator!=(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(b);
	return c.lower() != a.lower() || c.upper() != a.upper();
	}

// Binary != operator
// Works for all mixed classes
//
template<class _IT, class _X> inline bool operator!=(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);
	return c.lower() != b.lower() || c.upper() != b.upper();
	}

// != operator
// Works for all classes
//
template<class _IT, class _X> inline bool operator!=(const interval<_IT>& a, const interval<_X>& b)
	{
	return a.lower() != b.lower() || a.upper() != b.upper();
	}

// != operator
// Works for all classes
//
//template<class _IT> inline bool operator!=(const interval<_IT>& a, const interval<_IT>& b)
//	{
//	return a.lower() != b.lower() || a.upper() != b.upper();
//	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Boolean operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval abs()
///
//////////////////////////////////////////////////////////////////////////////////////

template<class _IT> inline interval<_IT> abs( const interval<_IT>& a )
	{
	if (a.lower() >= _IT(0) )
		return a;
	else
		if (a.upper() <= _IT(0) )
			return -a;

	return interval<_IT>(_IT(0), ( a.upper() > -a.lower() ? a.upper() : -a.lower() ) );
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END interval functions
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sqrt(), log10(), log(), exp() and pow()
///
//////////////////////////////////////////////////////////////////////////////////////

// Support function for correctly converting and double number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//

// Support function for correctly converting and double number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline float tofloat(const double& d, enum round_mode rm)
	{
	float fres;

	switch (rm)
		{
		case ROUND_DOWN:  fpdown();  break;
		case ROUND_UP: fpup();  break;
		}

	fres = (float)d;
	fpnear();
	return fres;
	}

// If H/W Support allows us to control the rounding mode then we can do it directly.
// Log(2) for double
inline double ln2double(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fldln2;				 Load ln2
		fstp qword ptr[res]; Store result in res
		}
	fpnear();

	return res;
	}

// Log(10) for double
inline double ln10double(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		; ln10 = FLDL2T * FLDLN2
		fldl2t;                 Load log2(10)
		fldln2;					Load LN2
		fmulp st(1),st;			Calculate LN2 * lOG2(10)
		fstp qword ptr[res];	Store ln10 in result
		}
	fpnear();

	return res;
	}

// PI for double
inline double pidouble(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fldpi;				 Load PI
		fstp qword ptr[res]; Store result in res
		}
	fpnear();

	return res;
	}

// Sqrt() for double
inline double sqrtdouble( double d, enum round_mode rm )
	{
	double sq = d;

	switch( rm )
	{
	case ROUND_DOWN: fpdown(); break;
	case ROUND_UP: fpup(); break;
	}

	_asm 
		{
		fld qword ptr[sq];  Load lower into floating point stack
		fsqrt;              Calculate sqrt
		fstp qword ptr[sq]; Store result in sq
		}
	fpnear();

	return sq;
	}

// Sqrt() for float
inline float sqrtfloat(float f, enum round_mode rm)
	{
	double sq = f;
	float fres;

	sq = sqrtdouble(sq, rm);
	fres = tofloat(sq, rm);
	return fres;
	}

// log() for double
inline double logdouble(double d, enum round_mode rm)
	{
	double lg = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[lg];      Load lower into floating point stack
		fldln2;                 Load loge2
		fxch st(1);             Exchange stack top
		fyl2x;                  Calculate y * ln2 x
		fstp qword ptr[lg];     Store result in lower
		}

	fpnear();

	return lg;
	}

// log() for float
inline float logfloat(float f, enum round_mode rm)
	{
	double lg = (double)f;
	float fres;

	lg = logdouble( lg, rm );
	fres = tofloat(lg, rm);
	return fres;
	}

// log10() for double
inline double log10double(double d, enum round_mode rm)
	{
	double lg = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[lg];      Load lower into floating point stack
		fldlg2;                 Load log10(2)
		fxch st(1);             Exchange stack top
		fyl2x;                  Calculate y * ln2 x
		fstp qword ptr[lg];     Store result in lg
		}

	fpnear();

	return lg;
	}

// log10 for float
inline float log10float(float f, enum round_mode rm)
	{
	double lg = (double)f;
	float fres;

	lg = log10double( lg, rm );
	fres = tofloat(lg, rm);
	return fres;
	}



// sqrt for float using managed code.
//
inline interval<float> sqrt( const interval<float>& x )
   {
   float lower, upper;
#ifdef HARDWARE_SUPPORT
   lower = sqrtfloat( x.lower(), ROUND_DOWN );
   upper = sqrtfloat( x.upper(), ROUND_UP );
#else
#endif
   return interval<float>( lower, upper );
   }

// sqrt for double using managed code.
//
inline interval<double> sqrt( const interval<double>& x )
   { 
   double lower, upper;

	lower = sqrtdouble( x.lower(), ROUND_DOWN );
	upper = sqrtdouble( x.upper(), ROUND_UP );
   return interval<double>( lower, upper );
   }


// log for float using manged code.
//
inline interval<float> log( const interval<float>& x )
	{
	float lower, upper;
 
	lower = logfloat( x.lower(), ROUND_DOWN);
	upper = logfloat( x.upper(), ROUND_UP);
   return interval<float>( lower, upper );
   }

// log for double using managed code.
//
inline interval<double> log( const interval<double>& x )
	{
	double lower, upper;

	lower = logdouble( x.lower(), ROUND_DOWN);
	upper = logdouble( x.upper(), ROUND_UP);
   return interval<double>( lower, upper );
   }

// log10 for float using manged code.
//
inline interval<float> log10( const interval<float>& x )
   {
   float lower, upper;

   lower = log10float( x.lower(), ROUND_DOWN);
   upper = log10float(x.upper(), ROUND_UP);
   return interval<float>( lower, upper );
   }

// log10 for double using managed code.
//
inline interval<double> log10( const interval<double>& x )
   {   
	double lower, upper;

	lower = log10double( x.lower(), ROUND_DOWN);
	upper = log10double( x.upper(), ROUND_UP);

   return interval<double>( lower, upper );
   }


// MSC exp() does not allow rounding control
// So we have to do it manually
// Use a taylor series until their is no more change in the result
// exp(x) == 1 + x + x^2/2!+x^3/3!+....
// Equivalent with the same standard C function call
// use argument reduction via exp(x)=(exp(x/2^k)2^k	
// And use Brent enhancement using the double formula:
// expm(x)=exp(x)-1 && expm(2x)=expm(x)(2+expm(x)) on the backend to preseve
// loss of significance digits
//
inline interval<double> exp(const interval<double>& x)
	{
	int  i, k = 0;
	interval<double> c, res, p0, old;
	const interval<double> c1(1), c2(2);

	c = x; 
	if (x.is_class() == NEGATIVE)
		c = abs(c);

	// Determine Reduction factor
	k = int((log(2) + log(abs(c.center()))) / log(2));
	k = std::min(k, 10);
	if (k > 0)
		{
		i = 2 << (k - 1);
		c /= interval<double>( (double)i );
		}

	p0 = c;
	old = c1;
	res = old + p0;
	for (i = 2; i < 100 && (res.lower() != old.lower() && res.upper() != old.upper()); i++)
		{
		old = res;
		p0 *= ( c / interval<double>((double)i));
		res += p0;
		}
	
	// Brent enhancement avoid loss of significant digits when x is small.
	if (k>0)
		{
		res -= c1;
		for (; k > 0; k--)
			res = (c2 + res)*res;
		res += c1;
		}

	if (x.is_class() == NEGATIVE)
		res = c1 / res;
	
	return res;
	}

// MSC exp() does not allow rounding control for the exp()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of exp() and convert back to float preserving as much accuracy as possible
//
inline interval<float> exp(const interval<float>& x)
	{
	interval<double> exp(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = exp(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<double>(lower, upper);
	}


// MSC pow() does not allow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
// To avoid loss of precision we actually perform the operation using double and then 
// convert the result back to float. This is consistent with the pow() that only takes double as an argument.
// 
inline interval<float> pow(const interval<float>& x, const float y)
	{
	interval<double> c(x);
	float upper, lower;

	c = log(c);
	c *= interval<double>(y);
	c = exp(c);
	lower = (float)c.lower();
	upper = (float)c.upper();
	if (lower > c.lower() )
		lower -= lower * 0.5f * FLT_EPSILON;
	if (upper < c.upper() )
		upper += upper * 0.5f * FLT_EPSILON;

	return interval<float>(lower, upper);
	}

// MSC pow() does not alllow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
// 
inline interval<double> pow(const interval<double>& x, const double y)
	{
	interval<double> c;

	c = log(x);
	c *= interval<double>(y);
	c = exp(c);

	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sqrt(), log10(), log(), exp(), pow()
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval constants like PI, LN2 and LN10
///
//////////////////////////////////////////////////////////////////////////////////////

// Load manifest constant PI for double
//
inline interval<double> int_pidouble()
	{
	interval<double> pi;
 
	pi.lower( pidouble(ROUND_DOWN) );
	pi.upper( pidouble(ROUND_UP) );
	return pi;
	}

// Load manifest constant PI for float
//
inline interval<float> int_pifloat()
	{
	interval<double> pid;
	interval<float> pif;

	pid = int_pidouble();
	pif.lower(tofloat(pid.lower(), ROUND_DOWN));
	pif.upper(tofloat(pid.upper(), ROUND_UP));
	return pif;;
	}



// Load manifest constant lN 2 for double
//
inline interval<double> int_ln2double()
	{
	interval<double> ln2;
 
	ln2.lower(ln2double(ROUND_DOWN));
	ln2.upper(ln2double(ROUND_UP));
	return ln2;
	}

// Load manifest constant LN2 for float
//
inline interval<float> int_ln2float()
	{
	interval<double> ln2d;
	interval<float> ln2f;

	ln2d = int_ln2double();
	ln2f.lower(tofloat(ln2d.lower(), ROUND_DOWN));
	ln2f.upper(tofloat(ln2d.upper(), ROUND_UP));
	return ln2f;;
	}

// Load manifest constant ln10 for double
//
inline interval<double> int_ln10double()
	{
	interval<double> ln10;
 
	ln10.lower(ln10double(ROUND_DOWN));
	ln10.upper(ln10double(ROUND_UP));
	return ln10;
	}

// Load manifest constant LN10 for float
//
inline interval<float> int_ln10float()
	{
	interval<double> ln10d;
	interval<float> ln10f;

	ln10d = int_ln10double();
	ln10f.lower(tofloat(ln10d.lower(), ROUND_DOWN));
	ln10f.upper(tofloat(ln10d.upper(), ROUND_UP));
	return ln10f;;
	}



//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval constants
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////


// Sin() for double interval
//
inline double sindouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fsin;				 Calculate sin
		fstp qword ptr[res]; Store result in lower
		}
	fpnear();

	return res;
	}

// Sin() for float interval
//
inline float sinfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = sindouble(res, rm);
	fres = tofloat(res, rm);
	return fres;
	}


// cos() for double interval
//
inline double cosdouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fcos;				 Calculate cos
		fstp qword ptr[res]; Store result in lower
		}
	fpnear();

	return res;
	}

// Cos(x) for float interval
//
inline float cosfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = cosdouble( res, rm );
	fres = tofloat(res, rm);
	return fres;
	}

// tan() for double interval
//
inline double tandouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fptan;               Calculate tan
		fstp qword ptr[res]; Pop ST(0) and ignore
		fstp qword ptr[res]; Store result
		}
	fpnear();

	return res;
	}

// Tan(x) for float interval
//
inline float tanfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = tandouble( res, rm );
	fres = tofloat( res, rm );
	return fres;
	}

// atan() for double interval
//
inline double atandouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];		Load lower into floating point stack
		fld1;					Load 1.0 on top of stack
		fpatan;					Calculate tan			
		fstp qword ptr[res];	Store result
		}	
	fpnear();

	return res;
	}

// Return atan(x) rounded up or rounded down
//
inline float atanfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = atandouble( res, rm );
	fres = tofloat(res, rm);
	return fres;  // or alternatively return tofloat(atandouble((double)f, rm ), rm );
	}


// Interval sin(x) for float
//
inline interval<float> sin(const interval<float>& x)
	{
	float lower, upper;

	lower = sinfloat(x.lower(), ROUND_DOWN);
	upper = sinfloat(x.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Interval Sin(x) for double
//
inline interval<double> sin(const interval<double>& x)
	{
	double lower, upper;

	lower = sindouble(x.lower(), ROUND_DOWN);
	upper = sindouble(x.upper(), ROUND_UP);
	return interval<double>(lower, upper);
	}



// Interval cos(x) for float.
//
inline interval<float> cos(const interval<float>& x)
	{
	float lower, upper;
 
	lower = cosfloat(x.lower(), ROUND_DOWN);
	upper = cosfloat(x.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Interval Cos(x) for double.
//
inline interval<double> cos(const interval<double>& x)
	{
	double lower, upper;
 
	lower = cosdouble(x.lower(), ROUND_DOWN);
	upper = cosdouble(x.upper(), ROUND_UP);
	return interval<double>(lower, upper);
	}



// Interval tan(x) for float.
//
inline interval<float> tan(const interval<float>& x)
	{
	float lower, upper;
 
	lower = tanfloat(x.lower(), ROUND_DOWN);
	upper = tanfloat(x.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Interval Tan(x) for double.
//
inline interval<double> tan(const interval<double>& x)
	{
	double lower, upper;
 
	lower = tandouble(x.lower(), ROUND_DOWN);
	upper = tandouble(x.upper(), ROUND_UP);
	return interval<double>(lower, upper);
	}

// Interval Arctan(x) for float.
//
inline interval<float> atan(const interval<float>& x)
	{
	float lower, upper;
 
	lower = atanfloat(x.lower(), ROUND_DOWN);
	upper = atanfloat(x.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}



// Interval ArcTan(x) for double.
//
inline interval<double> atan(const interval<double>& x)
	{
	double lower, upper;

	lower = atandouble(x.lower(), ROUND_DOWN);
	upper = atandouble(x.upper(), ROUND_UP);
	return interval<double>(lower, upper);
	}

// MSC asin() does not allow rounding control
// So we have to do it manually
/// Description:
///   Use a taylor series until their is no more change in the result
///   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
///   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	  This function replace the other function using newton iteration. Taylor series is significant
//    faster e.g 50% for 10 digits, 3x for 100 digits and 5x for 1000 digits.
//
inline interval<double> asin(const interval<double>& x)
	{
	int k, sign;
	interval<double> r, u, v, v2, sqrt2, lc, uc;
	const double c1(1), c2(2);

	if (x.lower() >= c1 || x.upper() <= -c1)
		{
		throw interval<double>::domain_error(); return x;
		}

	v = x;
	if (v.lower() < -c1)
		v.lower(-c1);
	if (v.upper() > c1)
		v.upper(c1);

	sign = v.is_class();
	if (sign == NEGATIVE)
		v = -v;

	// Now use the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
	// until argument is less than dlimit
	// Reduce the argument to below 0.5 to make the newton run faster
	sqrt2 = interval<double>(c2);				// Ensure correct number of digits
	sqrt2 = sqrt(sqrt2);	// Now calculate sqrt2 with precision digits
	for (k = 0; v.lower() > 0.5; k++)
		v /= sqrt2 * sqrt(interval<double>(c1) + sqrt(interval<double>(c1) - v * v));

	v2 = v * v;
	r = v;
	u = v;
	// Now iterate using taylor expansion
	for (unsigned int j = 3;; j += 2)
		{
		uc = interval<double>((j - 2) * (j - 2));
		lc = interval<double>(j * j - j);
		v = uc * v2 / lc;
		r *= v;
		if( u.lower() + r.lower() == u.lower() || u.upper() + r.upper() == u.upper())
			break;
		u += r;
		}

	if (k > 0)
		u *= interval<double>(1 << k);

	if (sign == NEGATIVE )
		u = -u;

	return u;
	}

// MSC acos() does not allow rounding control
// So we have to do it manually
/// Description:
///   Use a taylor series until their is no more change in the result
///   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
///   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	  This function replace the other function using newton iteration. Taylor series is significant
//    faster e.g 50% for 10 digits, 3x for 100 digits and 5x for 1000 digits.
//
inline interval<double> acos(const interval<double>& x)
	{
	interval<double> pi, res;
	const double c1(1);

	if (x.lower() >= c1 || x.upper() <= -c1)
		{
		throw interval<double>::domain_error(); return x;
		}

	pi.lower( pidouble(ROUND_DOWN) );
	pi.upper( pidouble(ROUND_UP) );
	res = pi * interval<double>(0.5) - asin( x );
	return res;
	}

// MSC asin() does not allow rounding control for the asin()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of asin() and convert back to float preserving as much accuracy as possible
//
inline interval<float> asin(const interval<float>& x)
	{
	interval<double> asin(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = asin(fx);
	lower = tofloat(fx.lower(), ROUND_DOWN);
	upper = tofloat(fx.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}


// MSC acos() does not allow rounding control for the acos()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of acos() and convert back to float preserving as much accuracy as possible
//
inline interval<float> acos(const interval<float>& x)
	{
	interval<double> acos(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = acos(fx);
	lower = tofloat(fx.lower(), ROUND_DOWN);
	upper = tofloat(fx.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}


//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////

#endif

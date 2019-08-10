
#ifndef _E_FLOAT_COMPLEX_CLR_2010_01_06_H_
  #define _E_FLOAT_COMPLEX_CLR_2010_01_06_H_

  #include <functions/complex/e_float_complex.h>
  #include <functions/elementary/elementary_complex.h>
  #include <interop/clr/e_float_template_clr.h>
  #include <interop/clr/e_float_clr.h>

  namespace e_float_clr
  {
    namespace cli
    {
      public ref class ef_complex sealed : public cli::e_float_template< ::ef_complex>
      {
      public:

        explicit cli::ef_complex(const  INT32 n) : e_float_template< ::ef_complex>(n) { }
        explicit cli::ef_complex(const  INT64 n) : e_float_template< ::ef_complex>(n) { }
        explicit cli::ef_complex(const UINT32 u) : e_float_template< ::ef_complex>(u) { }
        explicit cli::ef_complex(const UINT64 u) : e_float_template< ::ef_complex>(u) { }
        explicit cli::ef_complex(const double d) : e_float_template< ::ef_complex>(d) { }

        cli::ef_complex()                                   { delete my_p; my_p = new ::ef_complex(); }
        cli::ef_complex(cli::e_float^ re)                   { delete my_p; my_p = new ::ef_complex(re); }
        cli::ef_complex(cli::e_float^ re, cli::e_float^ im) { delete my_p; my_p = new ::ef_complex(re, im); }

        explicit cli::ef_complex(const ::ef_complex& f) : e_float_template< ::ef_complex>(f) { }
        cli::ef_complex(cli::ef_complex^ f) : e_float_template< ::ef_complex>(f) { }

        virtual cli::ef_complex::~ef_complex() { }

        System::String^ get_str(void) { return gcnew System::String((my_ref().get_str()).c_str()); }

        // Implement some of IronPython's overridden operators.
        cli::ef_complex^ __pos__  (void) { return gcnew cli::ef_complex(my_ref().get_pos()); }
        cli::ef_complex^ __neg__  (void) { return gcnew cli::ef_complex(my_ref().get_neg()); }
        cli::ef_complex^ __abs__  (void) { return gcnew cli::ef_complex(my_ref().get_abs()); }
        INT64            __int__  (void) { return my_ref().get_int(); }
        double           __float__(void) { return my_ref().get_flt(); }
        System::String^  __str__  (void) { return get_str(); }

        cli::e_float^ real(void) { return gcnew cli::e_float(my_ref().real()); }
        cli::e_float^ imag(void) { return gcnew cli::e_float(my_ref().imag()); }

        static cli::e_float^ real(cli::ef_complex^ z) { return z->real(); }
        static cli::e_float^ imag(cli::ef_complex^ z) { return z->imag(); }

        cli::e_float^ norm(void) { return gcnew cli::e_float(my_ref().norm()); }

        cli::ef_complex^ operator =(cli::e_float^ v) { const_cast< ::ef_complex&>(my_ref())  = (v->my_ref()); return this; }
        cli::ef_complex^ operator+=(cli::e_float^ v) { const_cast< ::ef_complex&>(my_ref()) += (v->my_ref()); return this; }
        cli::ef_complex^ operator-=(cli::e_float^ v) { const_cast< ::ef_complex&>(my_ref()) -= (v->my_ref()); return this; }
        cli::ef_complex^ operator*=(cli::e_float^ v) { const_cast< ::ef_complex&>(my_ref()) *= (v->my_ref()); return this; }
        cli::ef_complex^ operator/=(cli::e_float^ v) { const_cast< ::ef_complex&>(my_ref()) /= (v->my_ref()); return this; }

        cli::ef_complex^ operator+(void) { return this; }
        cli::ef_complex^ operator-(void) { return gcnew cli::ef_complex(-my_ref()); }

        // The CLR does not support postfix increment/decrement operators.
        // Implement explicit postfix increment/decrement functions.
        static cli::ef_complex^ post_inc(cli::ef_complex^ u) { cli::ef_complex^ v = gcnew cli::ef_complex(u); ++u; return v; }
        static cli::ef_complex^ post_dec(cli::ef_complex^ u) { cli::ef_complex^ v = gcnew cli::ef_complex(u); --u; return v; }

        static cli::ef_complex^ operator+(cli::ef_complex^ u, cli::ef_complex^ v) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) += v); }
        static cli::ef_complex^ operator-(cli::ef_complex^ u, cli::ef_complex^ v) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) -= v); }
        static cli::ef_complex^ operator*(cli::ef_complex^ u, cli::ef_complex^ v) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) *= v); }
        static cli::ef_complex^ operator/(cli::ef_complex^ u, cli::ef_complex^ v) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) /= v); }

        static cli::ef_complex^ operator+(cli::ef_complex^ u, cli::e_float^ v) { return gcnew cli::ef_complex(u) += v; }
        static cli::ef_complex^ operator-(cli::ef_complex^ u, cli::e_float^ v) { return gcnew cli::ef_complex(u) -= v; }
        static cli::ef_complex^ operator*(cli::ef_complex^ u, cli::e_float^ v) { return gcnew cli::ef_complex(u) *= v; }
        static cli::ef_complex^ operator/(cli::ef_complex^ u, cli::e_float^ v) { return gcnew cli::ef_complex(u) /= v; }

        static cli::ef_complex^ operator+(cli::e_float^ u, cli::ef_complex^ v) { return gcnew cli::ef_complex(::e_float(u) + (v->my_ref())); }
        static cli::ef_complex^ operator-(cli::e_float^ u, cli::ef_complex^ v) { return gcnew cli::ef_complex(::e_float(u) - (v->my_ref())); }
        static cli::ef_complex^ operator*(cli::e_float^ u, cli::ef_complex^ v) { return gcnew cli::ef_complex(::e_float(u) * (v->my_ref())); }
        static cli::ef_complex^ operator/(cli::e_float^ u, cli::ef_complex^ v) { return gcnew cli::ef_complex(::e_float(u) / (v->my_ref())); }

        static cli::ef_complex^ operator+(cli::ef_complex^ u, const INT32 n) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) += n); }
        static cli::ef_complex^ operator-(cli::ef_complex^ u, const INT32 n) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) -= n); }
        static cli::ef_complex^ operator*(cli::ef_complex^ u, const INT32 n) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) *= n); }
        static cli::ef_complex^ operator/(cli::ef_complex^ u, const INT32 n) { return static_cast<cli::ef_complex^>(gcnew cli::ef_complex(u) /= n); }

        static cli::ef_complex^ operator+(const INT32 n, cli::ef_complex^ u) { return gcnew cli::ef_complex(::e_float(n) + (u->my_ref())); }
        static cli::ef_complex^ operator-(const INT32 n, cli::ef_complex^ u) { return gcnew cli::ef_complex(::e_float(n) - (u->my_ref())); }
        static cli::ef_complex^ operator*(const INT32 n, cli::ef_complex^ u) { return gcnew cli::ef_complex(::e_float(n) * (u->my_ref())); }
        static cli::ef_complex^ operator/(const INT32 n, cli::ef_complex^ u) { return gcnew cli::ef_complex(::e_float(n) / (u->my_ref())); }

        static bool operator==(cli::ef_complex^ u, cli::ef_complex^ v) { return ((u->my_ref()) == (v->my_ref())); }
        static bool operator!=(cli::ef_complex^ u, cli::ef_complex^ v) { return ((u->my_ref()) != (v->my_ref())); }

        static bool operator==(cli::ef_complex^ u, cli::e_float^ v) { return ((u->my_ref()) == (v->my_ref())); }
        static bool operator!=(cli::ef_complex^ u, cli::e_float^ v) { return ((u->my_ref()) != (v->my_ref())); }

        static bool operator==(cli::e_float^ u, cli::ef_complex^ v) { return ((u->my_ref()) == (v->my_ref())); }
        static bool operator!=(cli::e_float^ u, cli::ef_complex^ v) { return ((u->my_ref()) != (v->my_ref())); }

        static bool operator==(cli::ef_complex^ u, const INT32 n) { return ((u->my_ref()) == n); }
        static bool operator!=(cli::ef_complex^ u, const INT32 n) { return ((u->my_ref()) != n); }

        static bool operator==(const INT32 n, cli::ef_complex^ u) { return (n == (u->my_ref())); }
        static bool operator!=(const INT32 n, cli::ef_complex^ u) { return (n != (u->my_ref())); }
      };
    }
  }

#endif // _E_FLOAT_COMPLEX_CLR_2010_01_06_H_

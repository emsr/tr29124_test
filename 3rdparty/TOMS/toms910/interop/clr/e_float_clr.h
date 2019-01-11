
#ifndef _E_FLOAT_CLR_2010_01_05_H_
  #define _E_FLOAT_CLR_2010_01_05_H_

  #include <e_float/e_float.h>
  #include <interop/clr/e_float_template_clr.h>

  namespace e_float_clr
  {
    namespace cli
    {
      public ref class e_float sealed : public cli::e_float_template< ::e_float>
      {
      public:

        explicit cli::e_float(const  INT32 n) : e_float_template< ::e_float>(n) { }
        explicit cli::e_float(const  INT64 n) : e_float_template< ::e_float>(n) { }
        explicit cli::e_float(const UINT32 u) : e_float_template< ::e_float>(u) { }
        explicit cli::e_float(const UINT64 u) : e_float_template< ::e_float>(u) { }
        explicit cli::e_float(const double d) : e_float_template< ::e_float>(d) { }

        explicit cli::e_float(const char* const s) { delete my_p; my_p = new ::e_float(s); }
        explicit cli::e_float(System::String^ sys_str);

        explicit cli::e_float(const ::e_float& f) : e_float_template< ::e_float>(f) { }
        cli::e_float(cli::e_float^ f) : e_float_template< ::e_float>(f) { }

        virtual cli::e_float::~e_float() { }

        System::String^ get_str(void) { return gcnew System::String((my_ref().get_str()).c_str()); }

        // Implement some of IronPython's overridden operators.
        cli::e_float^   __pos__  (void) { return gcnew cli::e_float(my_ref().get_pos()); }
        cli::e_float^   __neg__  (void) { return gcnew cli::e_float(my_ref().get_neg()); }
        cli::e_float^   __abs__  (void) { return gcnew cli::e_float(my_ref().get_abs()); }
        INT64           __int__  (void) { return my_ref().get_int(); }
        double          __float__(void) { return my_ref().get_flt(); }
        System::String^ __str__  (void) { return get_str(); }

        bool isnan   (void) { return my_ref().isnan   (); }
        bool isinf   (void) { return my_ref().isinf   (); }
        bool isfinite(void) { return my_ref().isfinite(); }
        bool iszero  (void) { return my_ref().iszero  (); }
        bool isone   (void) { return my_ref().isone   (); }
        bool isint   (void) { return my_ref().isint   (); }
        bool isneg   (void) { return my_ref().isneg   (); }

        INT64 order(void) { return my_ref().order(); }

        cli::e_float^ operator+(void) { return this; }
        cli::e_float^ operator-(void) { return gcnew cli::e_float(-my_ref()); }

        // The CLR does not support postfix increment/decrement operators.
        // Implement explicit postfix increment/decrement functions.
        static cli::e_float^ post_inc(cli::e_float^ u) { cli::e_float^ v = gcnew cli::e_float(u); ++u; return v; }
        static cli::e_float^ post_dec(cli::e_float^ u) { cli::e_float^ v = gcnew cli::e_float(u); --u; return v; }

        static cli::e_float^ operator+(cli::e_float^ u, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) += v); }
        static cli::e_float^ operator-(cli::e_float^ u, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) -= v); }
        static cli::e_float^ operator*(cli::e_float^ u, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) *= v); }
        static cli::e_float^ operator/(cli::e_float^ u, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) /= v); }

        static cli::e_float^ operator+(cli::e_float^ u, const INT32 n) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) += n); }
        static cli::e_float^ operator-(cli::e_float^ u, const INT32 n) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) -= n); }
        static cli::e_float^ operator*(cli::e_float^ u, const INT32 n) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) *= n); }
        static cli::e_float^ operator/(cli::e_float^ u, const INT32 n) { return static_cast<cli::e_float^>(gcnew cli::e_float(u) /= n); }

        static cli::e_float^ operator+(const INT32 n, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(n) += v); }
        static cli::e_float^ operator-(const INT32 n, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(n) -= v); }
        static cli::e_float^ operator*(const INT32 n, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(n) *= v); }
        static cli::e_float^ operator/(const INT32 n, cli::e_float^ v) { return static_cast<cli::e_float^>(gcnew cli::e_float(n) /= v); }

        static bool operator< (cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) <  (v->my_ref())); }
        static bool operator<=(cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) <= (v->my_ref())); }
        static bool operator==(cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) == (v->my_ref())); }
        static bool operator!=(cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) != (v->my_ref())); }
        static bool operator>=(cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) >= (v->my_ref())); }
        static bool operator> (cli::e_float^ u, cli::e_float^ v) { return ((u->my_ref()) >  (v->my_ref())); }
                                                                                                         
        static bool operator< (cli::e_float^ u, const INT32 n) { return ((u->my_ref()) <  n); }
        static bool operator<=(cli::e_float^ u, const INT32 n) { return ((u->my_ref()) <= n); }
        static bool operator==(cli::e_float^ u, const INT32 n) { return ((u->my_ref()) == n); }
        static bool operator!=(cli::e_float^ u, const INT32 n) { return ((u->my_ref()) != n); }
        static bool operator>=(cli::e_float^ u, const INT32 n) { return ((u->my_ref()) >= n); }
        static bool operator> (cli::e_float^ u, const INT32 n) { return ((u->my_ref()) >  n); }
                                                                                           
        static bool operator< (const INT32 n, cli::e_float^ u) { return (n <  (u->my_ref())); }
        static bool operator<=(const INT32 n, cli::e_float^ u) { return (n <= (u->my_ref())); }
        static bool operator==(const INT32 n, cli::e_float^ u) { return (n == (u->my_ref())); }
        static bool operator!=(const INT32 n, cli::e_float^ u) { return (n != (u->my_ref())); }
        static bool operator>=(const INT32 n, cli::e_float^ u) { return (n >= (u->my_ref())); }
        static bool operator> (const INT32 n, cli::e_float^ u) { return (n >  (u->my_ref())); }
      };                                                                                   
    }
  }

#endif // _E_FLOAT_CLR_2010_01_05_H_

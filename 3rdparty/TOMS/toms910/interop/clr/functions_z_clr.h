
#ifndef _FUNCTIONS_Z_CLR_2010_01_18_H_
  #define _FUNCTIONS_Z_CLR_2010_01_18_H_

  #include <functions/functions.h>
  #include <interop/clr/e_float_clr.h>
  #include <interop/clr/e_float_complex_clr.h>
  #include <interop/clr/functions_base_clr.h>

  namespace e_float_clr
  {
    namespace cli
    {
      public value struct efz sealed : public functions_base
      {
        // <functions/elementary/elementary_complex.h>
        static cli::e_float^ norm(cli::ef_complex^ z) { return z->norm(); }
        static cli::e_float^ abs (cli::ef_complex^ z) { return gcnew cli::e_float(::efz::abs(z)); }
        static cli::e_float^ arg (cli::ef_complex^ z) { return gcnew cli::e_float(::efz::arg(z)); }
        static cli::e_float^ real(cli::ef_complex^ z) { return z->real(); }
        static cli::e_float^ imag(cli::ef_complex^ z) { return z->imag(); }

        static cli::ef_complex^ conj(cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::conj(z)); }
        static cli::ef_complex^ iz  (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::iz(z)); }

        static cli::ef_complex^ polar(cli::e_float^ mod, cli::e_float^ arg) { return gcnew cli::ef_complex(::efz::polar(mod, arg)); }

        static cli::ef_complex^ sin     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::sin(z)); }
        static cli::ef_complex^ cos     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::cos(z)); }
        static cli::ef_complex^ tan     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::tan(z)); }
        static void             sincos  (cli::ef_complex^ z, cli::ef_complex^ s, cli::ef_complex^ c) { ::efz::sincos(z, const_cast< ::ef_complex*>(s->my_ptr()), const_cast< ::ef_complex*>(c->my_ptr())); }
        static cli::ef_complex^ csc     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::csc  (z)); }
        static cli::ef_complex^ sec     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::sec  (z)); }
        static cli::ef_complex^ cot     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::cot  (z)); }
        static cli::ef_complex^ asin    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::asin (z)); }
        static cli::ef_complex^ acos    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::acos (z)); }
        static cli::ef_complex^ atan    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::atan (z)); }
        static cli::ef_complex^ inv     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::inv  (z)); }
        static cli::ef_complex^ sqrt    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::sqrt (z)); }
        static cli::ef_complex^ exp     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::exp  (z)); }
        static cli::ef_complex^ log     (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::log  (z)); }
        static cli::ef_complex^ log10   (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::log10(z)); }
        static cli::ef_complex^ loga    (cli::ef_complex^ a, cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::loga(a, z)); }
        static cli::ef_complex^ pown    (cli::ef_complex^ z, const INT64 p)      { return gcnew cli::ef_complex(::efz::pown(z, p)); }
        static cli::ef_complex^ pow     (cli::ef_complex^ z, cli::ef_complex^ a) { return gcnew cli::ef_complex(::efz::pow(z, a)); }
        static cli::ef_complex^ rootn   (cli::ef_complex^ z, const INT32 p)      { return gcnew cli::ef_complex(::efz::rootn(z, p)); }
        static cli::ef_complex^ sinh    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::sinh(z)); }
        static cli::ef_complex^ cosh    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::cosh(z)); }
        static cli::ef_complex^ tanh    (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::tanh(z)); }
        static void             sinhcosh(cli::ef_complex^ z, cli::ef_complex^ s, cli::ef_complex^ c) { ::efz::sinhcosh(z, const_cast< ::ef_complex*>(s->my_ptr()), const_cast< ::ef_complex*>(c->my_ptr())); }
        static cli::ef_complex^ asinh   (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::asinh(z)); }
        static cli::ef_complex^ acosh   (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::acosh(z)); }
        static cli::ef_complex^ atanh   (cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::atanh(z)); }

        // <functions/gamma/gamma.h>
        static cli::ef_complex^ gamma     (cli::ef_complex^ z)                     { return gcnew cli::ef_complex(::efz::gamma(z)); }
        static cli::ef_complex^ beta      (cli::ef_complex^ a, cli::ef_complex^ b) { return gcnew cli::ef_complex(::efz::beta(a, b)); }
        static cli::ef_complex^ pochhammer(cli::ef_complex^ z, const UINT32 n)     { return gcnew cli::ef_complex(::efz::pochhammer(z, n)); }
        static cli::ef_complex^ pochhammer(cli::ef_complex^ z, cli::ef_complex^ a) { return gcnew cli::ef_complex(::efz::pochhammer(z, a)); }

        // <functions/polynomials/polynomials.h>
        static cli::ef_complex^ chebyshev_t(const INT32 n, cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::chebyshev_t(n, z)); }
        static cli::ef_complex^ chebyshev_u(const INT32 n, cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::chebyshev_u(n, z)); }
        static cli::ef_complex^ hermite    (const INT32 n, cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::hermite    (n, z)); }
        static cli::ef_complex^ laguerre   (const INT32 n, cli::ef_complex^ z) { return gcnew cli::ef_complex(::efz::laguerre   (n, z)); }

        static cli::ef_complex^ chebyshev_t(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp);
        static cli::ef_complex^ chebyshev_u(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp);
        static cli::ef_complex^ hermite    (const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp);
        static cli::ef_complex^ laguerre   (const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp);

        // <functions/zeta/zeta.h>
        static cli::ef_complex^ riemann_zeta(cli::ef_complex^ s)                     { return gcnew cli::ef_complex(::efz::riemann_zeta(s)); }
        static cli::ef_complex^ hurwitz_zeta(cli::ef_complex^ s, const INT32 n)      { return gcnew cli::ef_complex(::efz::hurwitz_zeta(s, n)); }
        static cli::ef_complex^ hurwitz_zeta(cli::ef_complex^ s, cli::ef_complex^ a) { return gcnew cli::ef_complex(::efz::hurwitz_zeta(s, a)); }
      };
    }
  }

#endif // _FUNCTIONS_Z_CLR_2010_01_18_H_

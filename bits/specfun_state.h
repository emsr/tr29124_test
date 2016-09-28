
  template<typename _Tp>
    struct __airy_t
    {
      _Tp x;
      _Tp Ai;
      _Tp Aip;
      _Tp Bi;
      _Tp Bip;
    };

  template<typename _Tp>
    struct __cyl_bessel_t
    {
      _Tp x;
      _Tp nu;
      _Tp J;
      _Tp Jp;
      _Tp N;
      _Tp Np;
    };

  template<typename _Tp>
    struct __cyl_mod_bessel_t
    {
      _Tp x;
      _Tp nu;
      _Tp I;
      _Tp Ip;
      _Tp K;
      _Tp Kp;
    };

  template<typename _Tp>
    struct __cyl_hankel_t
    {
      _Tp x;
      _Tp nu;
      _Tp H1;
      _Tp H1p;
      _Tp H2;
      _Tp H2p;
    };

  template<typename _Tp>
    struct __sph_bessel_t
    {
      _Tp x;
      unsigned int n;
      _Tp j;
      _Tp jp;
      _Tp n;
      _Tp np;
    };

  template<typename _Tp>
    struct __sph_mod_bessel_t
    {
      _Tp x;
      unsigned int n;
      _Tp i;
      _Tp ip;
      _Tp k;
      _Tp kp;
    };

  template<typename _Tp>
    struct __sph_hankel_t
    {
      _Tp x;
      unsigned int n;
      _Tp h1;
      _Tp h1p;
      _Tp h2;
      _Tp h2p;
    };

  template<typename _Tp>
    struct __pqgamma_t
    {
    };

  template<typename _Tp>
    struct __lgamma_t
    {
      _Tp lgamma;
      int sign;
    };

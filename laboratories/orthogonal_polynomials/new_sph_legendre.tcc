
  template<typename Tp>
    Tp
    sph_legendre(const unsigned int l, const unsigned int m,
                   const Tp theta)
    {
      constexpr auto _S_pi = emsr::math_constants<Tp>::pi;
      const auto x = std::cos(theta);
      if (m < 0 || l < m || x < -Tp{1} || x > Tp{1})
        throw std::domain_error("sph_legendre: bad argument");
      else if (m == 0)
        {
          Tp _P_l = legendre_p(l, x);
          Tp fact = std::sqrt(Tp(2 * l + 1) / (Tp{4} * _S_pi));
          _P_l *= fact;
          return _P_l;
        }
      else if (x == Tp{1} || x == -Tp{1})
        return Tp{0};
      else
        {
          // Starting value for recursion.
          // Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) )
	  //            (-1)^m (1-x^2)^(m/2) / pi^(1/4)
          const auto sgn = Tp(m % 2 == 1 ? -Tp{1} : Tp{1});
          const auto _Y_mp1m_factor = x * std::sqrt(Tp(2 * m + 3));
          const auto lncirc = std::log1p(-x * x);
          // Gamma(m+1/2) / Gamma(m)
          const auto lnpoch = std::lgamma(Tp(m + 0.5L))
                              - std::lgamma(Tp(m));
          const auto lnpre_val = -Tp{0.25L} * _TR1_M_LNPI
				 + Tp{0.5L} * (lnpoch + m * lncirc);
          const auto sr = std::sqrt((Tp{2} + Tp{1} / m)
				    / (Tp{4} * _S_pi));
          auto _Y_mm = sgn * sr * std::exp(lnpre_val);
          auto _Y_mp1m = _Y_mp1m_factor * _Y_mm;

          if (l == m)
            return _Y_mm;
          else if (l == m + 1)
            return _Y_mp1m;
          else
            {
              Tp _Y_lm = Tp{0};

              // Compute Y_l^m, l > m+1, upward recursion on l.
              for (int ell = m + 2; ell <= l; ++ell)
                {
                  const auto rat1 = Tp(ell - m) / Tp(ell + m);
                  const auto rat2 = Tp(ell - m - 1) / Tp(ell + m - 1);
                  const auto factor1 = std::sqrt(rat1
						 * Tp(2 * ell + 1)
						 * Tp(2 * ell - 1));
                  const auto factor2 = std::sqrt(rat1 * rat2
						 * Tp(2 * ell + 1)
						 / Tp(2 * ell - 3));
                  _Y_lm = (x * _Y_mp1m * factor1
			- Tp(ell + m - 1) * _Y_mm * factor2)
			/ Tp(ell - m);
                  _Y_mm = _Y_mp1m;
                  _Y_mp1m = _Y_lm;
                }

              return _Y_lm;
            }
        }
    }

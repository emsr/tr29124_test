

  /**
   *  @brief Evaluate polynomial based on the confluent hypergeometric
   *         representation.
   *  @f[
   *    L_n^\alpha(x) = (\alpha + 1)_n / n! _1F_1(-n, \alpha + 1, x)
   *  @f]
   *  Assumes n > 0 and @f$ \alpha @f$ != negative integer greater than -n.
   */
  template<typename _Tpa, class Tp>
  Tp
    laguerre_n_cp(unsigned int n, _Tpa alpha, Tp x)
    {
      constexpr auto max = Tp{0.9L} * std::numeric_limits<Tp>::max();
      const auto lnfact = log_gamma(n + Tp{1});
      const auto lg1 = log_gamma(alpha + Tp{1} + n);
      const auto s1 = log_gamma_sign(alpha + Tp{1} + n);
      const auto lg2 = log_gamma(alpha + Tp{1});
      const auto s2 = log_gamma_sign(alpha + Tp{1});
      const auto lnpre = (lg1 - lg2) - lnfact;

      auto poly_1F1 = Tp{1};
      for (auto k = int(n) - 1; k >= 0; --k)
	{
	  auto t = (-int(n) + k) * (x / Tp(k + 1))
		   / (alpha + Tp{1} + k);
	  auto r = t + Tp{1} / poly_1F1;
	  if (r > max / poly_1F1) // Internal error only, catch in the main routine.
            throw std::runtime_error("laguerre_n_cp: series failed);
	  else // Collect the Horner terms.
            poly_1F1  = Tp{1} + t * poly_1F1;
	}

      return std::exp(lnpre) * poly_1F1;
    }

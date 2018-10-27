// 
// Computes the hypergeometric function \mathbf{M}(a;c;z), using a Taylor
// series (Method (b)), up to a tolerance tol.
// 
template<typename Tp>
  Tp
  __hypfun_M_taylor_c(Tp a, Tp c, Tp z, Tp tol)
  {
    // Initialise vector r, which stores multipliers
    auto r = a / c;
    auto Akm2 = 1 + z * r;
    r = (a+1) / 2 / (c+1);
    Akm1 = Akm2 + z*z * a / c * r;

    for (int j = 3; j <= 500; ++j)
    {
      r = (a+j-1) / j / (c+j-1);
      // Compute A_j recursively
      Ak = Akm1 + (Akm1 - Akm2) * r * z;
      if (std::abs(Ak - Akm1) / std::abs(Akm1) < tol
       && std::abs(Akm1 - Akm2) / std::abs(Akm2) < tol)
        return Ak / std::tgamma(c);
      if (j == 500)
	std::__throw_runtime_error("__hypfun_M_taylor_c: no convergence");
    }
  }

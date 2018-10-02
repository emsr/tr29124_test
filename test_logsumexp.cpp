
template<class T>
  constexpr const T&
  logsumexp(const T& a, const T& b)
  {
    constexpr const auto& m = std::minmax(a, b);
    return m.second + std::log(T{1} + std::exp(m.first));
  }

template<class T, class Compare>
  constexpr const T&
  logsumexp(const T& a, const T& b, Compare comp)
  {
    constexpr const auto& m = std::minmax(a, b, comp);
    return m.second + std::log(T{1} + std::exp(m.first));
  }

template<class T>
  constexpr T
  logsumexp(std::initializer_list<T> ilist)
  {
    constexpr auto m = std::max(ilist);
    auto sum = T{};
    for (auto&& x : ilist)
      sum += std::exp(x - m);
    return m + std::log(sum);
  }

template<class T, class Compare>
  constexpr T
  logsumexp(std::initializer_list<T> ilist, Compare comp)
  {
    constexpr auto m = std::max(ilist, comp);
    auto sum = T{};
    for (auto&& x : ilist)
      sum += std::exp(x - m);
    return m + std::log(sum);
  }

// Similar functions...

// Safe sigmoid.
template<class T>
  constexpr T
  sigmoid(T x)
  {
    if (x >= T{0})
      return T{1} / (T{1} + std::exp(-x));
    else
      {
	constexpr const auto z = std::exp(x);
	return z / (T{1} + z);
      }
  }

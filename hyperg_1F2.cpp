
template<typename _Tp>
  _Tp
  __hyperg_1f2(_Tp a, _Tp b, _Tp c, _Tp x, _Tp& err)
  {
    auto an = a;
    auto bn = b;
    auto cn = c;
    auto a0 = _Tp{1};
    auto sum = _Tp{1};
    auto n = 1;
    auto t = _Tp{1};
    auto max = _Tp{0};

    do
      {
	if (an == _Tp{0})
	  break;
	if (bn == _Tp{0})
	  throw std::runtime_error("__hyperg_1f2: series failed");
	if (cn == _Tp{0})
	  throw std::runtime_error("__hyperg_1f2: series failed");
	if (a0 > 1.0e34 || n > 200)
	  throw std::runtime_error("__hyperg_1f2: series failed");
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += _Tp{1};
	bn += _Tp{1};
	cn += _Tp{1};
	++n;
	z = fabs(a0);
	if (z > max)
	  max = z;
	if (sum != _Tp{0})
	  t = fabs(a0 / sum);
	else
	  t = z;
      }
    while (t > stop);

    err = fabs(MACHEP * max / sum);

    return sum;
  }

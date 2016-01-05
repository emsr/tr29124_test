
template<typename _Tp>
  _Tp
  __hyperg_3f0(_Tp a, _Tp b, _Tp c, _Tp x, _Tp& err)
  {
    auto an = a;
    auto bn = b;
    auto cn = c;
    auto a0 = _Tp{1};
    auto sum = _Tp{1};
    auto n = 1;
    auto t = _Tp{1};
    auto max = _Tp{0};
    auto conv = _Tp{1.0e38};
    auto conv1 = conv;

    do
      {
	if (an == _Tp{0})
	  break;
	if (bn == _Tp{0})
	  break;
	if (cn == _Tp{0})
	  break;
	if (a0 > 1.0e34 || n > 200)
	  throw std::runtime_error("__hyperg_3f0: series failed");
	a0 *= (an * bn * cn * x) / n;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	++n;
	z = fabs(a0);
	if (z > max)
	  max = z;
	if (z >= conv)
	  {
	    if (z < max && z > conv1)
	      break;
	  }
	conv1 = conv;
	conv = z;
	sum += a0;
	if (sum != 0)
	  t = fabs(a0 / sum);
	else
	  t = z;
      }
    while (t > stop);

    t = fabs(MACHEP * max / sum);

    max = fabs(conv / sum);
    if (max > t)
      t = max;

    err = t;
    return sum;
  }

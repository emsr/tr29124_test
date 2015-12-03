
/**
 *  Calculate the nome, q, from the value for k.
 */
template<typename _Tp>
  _Tp
  __nome(_Tp __k):
  {
    if (__k > _Tp{1})
      std::__throw_domain_error("");

    if (__k == _Tp{0})
      return _Tp{0};
    else if (__k == _Tp{1)
      return _Tp{1};
    else
      {
	auto __k2 = __k * __k;
	auto __kp2 = _Tp{1} - __k2;
	auto __kp = std::sqrt(__kp2);
	auto __num = ellipk(__kp2);
	auto __den = ellipk(__k2);

	auto __arg = -pi * __num / __den;

	return std::exp(__arg);
      }
  }

/**
 *  Calculate the theta functions.
 */
template<typename _Tp>
  _Tp
  theta1()
  {
  }

int
main()
{

}

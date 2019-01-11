#include <TestCases/TestCaseDescription.h>

const TestCaseDescription& TestCaseDescription_case_00011_various_elem_math(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::floor(ef::pi()));"),
    std::string("        data.push_back(ef::ceil (ef::pi()));"),
    std::string("        data.push_back(ef::floor(-100 - ef::euler_gamma()));"),
    std::string("        data.push_back(ef::ceil (-100 - ef::euler_gamma()));"),
    std::string("        data.push_back(e_float(ef::to_int32(e_float(\"1e9\"))));"),
    std::string("        data.push_back(e_float(ef::to_int64(e_float(\"1e18\"))));"),
    std::string("        data.push_back(e_float(ef::to_int32(e_float(\"1e29\"))));"),
    std::string("        data.push_back(e_float(ef::to_int64(e_float(\"1e29\"))));"),
    std::string("        data.push_back(e_float(ef::to_int32(ef_complex(ef::pi(), ef::euler_gamma()))));"),
    std::string("        data.push_back(e_float(ef::to_int64(ef_complex(ef::pi(), ef::euler_gamma()))));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00011_various_elem_math",
                                       "{3, 4, -101, -100, 10^9, 10^18, 2147483647, 9223372036854775807, 3, 3}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00021_bernoulli(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::bernoulli(10u * k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00021_bernoulli",
                                       "Table[N[BernoulliB[10 k], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00022_euler(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::euler(10u * k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00022_euler",
                                       "Table[N[EulerE[10 k], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00023_stirling2(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const UINT32 nn = k;"),
    std::string("          const UINT32 kn = static_cast<UINT32>(k / static_cast<UINT32>(2u));"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::stirling2(nn, kn);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00023_stirling2",
                                       "Table[N[StirlingS2[k, IntegerPart[k / 2]], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00051_factorial(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::factorial(10u * k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00051_factorial",
                                       "Table[N[Factorial[10 k], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00052_factorial2(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::factorial2(10u * k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00052_factorial2",
                                       "Table[N[Factorial2[10 k], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00071_various_int_func(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::factorial2(-1301));"),
    std::string("        data.push_back(ef::factorial2(-2501));"),
    std::string("        data.push_back(ef::bernoulli(1u));"),
    std::string("        data.push_back(ef::bernoulli(10000u));"),
    std::string("        data.push_back(ef::bernoulli(10001u));"),
    std::string("        data.push_back(ef::bernoulli(10002u));"),
    std::string("        data.push_back(ef::euler(10000u));"),
    std::string("        data.push_back(ef::euler(10001u));"),
    std::string("        data.push_back(ef::euler(10002u));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00071_various_int_func",
                                       "{N[Factorial2[-1301], 400], N[Factorial2[-2501], 400], N[BernoulliB[1], 400], N[BernoulliB[10000], 400], N[BernoulliB[10001], 400], N[BernoulliB[10002], 400], N[EulerE[10000], 400], N[EulerE[10001], 400], N[EulerE[10002], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00101_sin(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::sin(ef::euler_gamma() * ((100 * k) - 5000));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00101_sin",
                                       "Table[N[Sin[EulerGamma ((100 k) - 5000)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00102_cos(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cos(ef::euler_gamma() * ((100 * k) - 5000));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00102_cos",
                                       "Table[N[Cos[EulerGamma ((100 k) - 5000)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00103_exp(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::exp(ef::sqrt((ef::pi() * (100 * k)) * (100 * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00103_exp",
                                       "Table[N[Exp[Sqrt[Pi (100 k) (100 k)]], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00104_log(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::log(ef::tenth() + ((ef::pi() * (100 * k)) * (100 * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00104_log",
                                       "Table[N[Log[(1/10) + (Pi (100 k) (100 k))], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00105_sqrt(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::sqrt(ef::pi() * k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00105_sqrt",
                                       "Table[N[Sqrt[Pi k], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00106_rootn(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k_plus_one = static_cast<INT32>(k + 1);"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::rootn(ef::googol() * (k_plus_one * ef::pi()), k_plus_one);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00106_rootn",
                                       "Table[N[Power[(10^100) (Pi (k + 1)), 1/(k + 1)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00111_sin_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::sin((ef::euler_gamma() + k) / 53);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00111_sin_small_x",
                                       "Table[N[Sin[(EulerGamma + k) / 53], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00112_cos_x_near_pi_half(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cos(ef::pi_half() - ((ef::euler_gamma() + k) / 523));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00112_cos_x_near_pi_half",
                                       "Table[N[Cos[(Pi / 2) - ((EulerGamma + k) / 523)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00113_atan_x_small_to_large(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_4k = e_float(\"1e-100\");"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::atan(ten_pow_4k);"),
    std::string("          ten_pow_4k *= static_cast<INT32>(10000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00113_atan_x_small_to_large",
                                       "Table[N[ArcTan[10^(-100 + (4 k))], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00114_various_trig(void)
{
  static const std::tr1::array<std::string, 3u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::sec(ef::euler_gamma()));"),
    std::string("        data.push_back(ef::sec(ef::catalan()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00114_various_trig",
                                       "{N[Sec[EulerGamma], 400], N[Sec[Catalan], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00115_various_elem_trans(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::sqrt(e_float(11.1)));"),
    std::string("        data.push_back(ef::sqrt(e_float(22.2)));"),
    std::string("        data.push_back(ef::log1p(static_cast<INT32>(-1) / ef::hundred()));"),
    std::string("        data.push_back(ef::log1p(static_cast<INT32>(+1) / ef::hundred()));"),
    std::string("        data.push_back(ef::log10(ef::googol() * ef::golden_ratio()));"),
    std::string("        data.push_back(ef::loga (ef::pi(), ef::googol() * ef::golden_ratio()));"),
    std::string("        data.push_back(ef::log  (e_float(\"1e123456\")));"),
    std::string("        data.push_back(ef::rootn(ef::euler_gamma(), static_cast<INT32>(-7)));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00115_various_elem_trans",
                                       "{N[Sqrt[111/10],400], N[Sqrt[222/10],400], N[Log[1 - (1 / 100)],400], N[Log[1 + (1 / 100)],400], N[Log[10, (10^100) GoldenRatio],400], N[Log[Pi, (10^100) GoldenRatio],400], N[Log[10^123456],400], N[Power[EulerGamma, -(1/7)], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00121_sinh(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::sinh(x * x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00121_sinh",
                                       "Table[N[Sinh[(EulerGamma + k)^2], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00122_cosh(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cosh(x * x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00122_cosh",
                                       "Table[N[Cosh[(EulerGamma + k)^2], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00123_tanh(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::tanh(x * x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00123_tanh",
                                       "Table[N[Tanh[(EulerGamma + k)^2], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00124_asinh(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::asinh(static_cast<INT32>(-25 + k) + ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00124_asinh",
                                       "Table[N[ArcSinh[(-25 + k) + EulerGamma], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00125_acosh(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::acosh(static_cast<INT32>(1 + k) + ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00125_acosh",
                                       "Table[N[ArcCosh[(1 + k) + EulerGamma], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00126_atanh(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::atanh((e_float(-25 + k) / static_cast<INT32>(26)) + (ef::euler_gamma() / static_cast<INT32>(100)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00126_atanh",
                                       "Table[N[ArcTanh[((-25 + k) / 26) + (EulerGamma / 100)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00201_gamma(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::gamma(ef::tenth() + (((10 * k) - 500) * ef::pi()));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00201_gamma",
                                       "Table[N[Gamma[(1/10) + (((10 k) - 500) Pi)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00202_gamma_medium_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::gamma(ef::third() + k);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00202_gamma_medium_x",
                                       "Table[N[Gamma[(1/3) + k], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00203_gamma_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        static const e_float fifty_three = e_float(53);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::gamma(ef::third() + (k / fifty_three));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00203_gamma_small_x",
                                       "Table[N[Gamma[(1/3) + (k / 53)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00204_gamma_tiny_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_two_k = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::gamma(ef::third() / ten_pow_two_k);"),
    std::string("          ten_pow_two_k *= 100;"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00204_gamma_tiny_x",
                                       "Table[N[Gamma[(1/3) (10^(-2 k))], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00205_gamma_near_neg_n(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_k = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::gamma_near_n(-87, (ef::third() / ten_pow_k));"),
    std::string("          ten_pow_k *= 10;"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00205_gamma_near_neg_n",
                                       "Table[N[Gamma[-87 + ((1/3) (10^(-k)))], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00221_various_gamma_func(void)
{
  static const std::tr1::array<std::string, 22u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::binomial(  20u,   10u));"),
    std::string("        data.push_back(ef::binomial( 200u,  100u));"),
    std::string("        data.push_back(ef::binomial(2000u, 1000u));"),
    std::string("        data.push_back(ef::binomial(  20u, ef::pi()));"),
    std::string("        data.push_back(ef::binomial( 200u, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::binomial(ef::pi(), 20u));"),
    std::string("        data.push_back(ef::binomial(ef::euler_gamma(), 200u));"),
    std::string("        data.push_back(ef::binomial(ef::pi(), ef::euler_gamma()));"),
    std::string("        data.push_back(efz::real(efz::pochhammer(ef_complex(ef::third(), ef::catalan()), ef_complex(ef::pi(), ef::euler_gamma()))));"),
    std::string("        data.push_back(efz::imag(efz::pochhammer(ef_complex(ef::third(), ef::catalan()), ef_complex(ef::pi(), ef::euler_gamma()))));"),
    std::string("        data.push_back(efz::real(efz::pochhammer(ef_complex(ef::third(), ef::catalan()), 17u)));"),
    std::string("        data.push_back(efz::imag(efz::pochhammer(ef_complex(ef::third(), ef::catalan()), 17u)));"),
    std::string("        data.push_back(ef::beta(20 + ef::pi(), 30 + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::incomplete_beta(ef::third(), 20 + ef::pi(), 30 + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::incomplete_beta(e_float(7) / 10, 13 + ef::pi(), 17 + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::incomplete_beta(e_float(99) / 100, 13 + ef::pi(), 17 + ef::euler_gamma()));"),
    std::string("        data.push_back(efz::real(efz::beta(ef_complex(ef::third(), ef::catalan()), ef_complex(ef::pi(), ef::euler_gamma()))));"),
    std::string("        data.push_back(efz::imag(efz::beta(ef_complex(ef::third(), ef::catalan()), ef_complex(ef::pi(), ef::euler_gamma()))));"),
    std::string("        data.push_back(ef::incomplete_gamma(ef::quarter() + ef::pi(), ef::third()));"),
    std::string("        data.push_back(ef::incomplete_gamma(101 + ef::pi(), e_float(31) / 3));"),
    std::string("        data.push_back(ef::gen_incomplete_gamma(ef::quarter() + ef::pi(), e_float(11) / 3, e_float(29) / 3));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00221_various_gamma_func",
                                       "{N[Binomial[20, 10], 400], N[Binomial[200, 100], 400], N[Binomial[2000, 1000], 400], N[Binomial[20, Pi], 400], N[Binomial[200, EulerGamma], 400], N[Binomial[Pi, 20], 400], N[Binomial[EulerGamma, 200], 400], N[Binomial[Pi, EulerGamma], 400], N[Re[Pochhammer[(1/3) + (Catalan I), Pi + (EulerGamma I)]], 400], N[Im[Pochhammer[(1/3) + (Catalan I), Pi + (EulerGamma I)]], 400], N[Re[Pochhammer[(1/3) + (Catalan I), 17]], 400], N[Im[Pochhammer[(1/3) + (Catalan I), 17]], 400], N[Beta[20 + Pi, 30 + EulerGamma], 400], N[Beta[1/3, 20 + Pi, 30 + EulerGamma], 400], N[Beta[7/10, 13 + Pi, 17 + EulerGamma], 400], N[Beta[99/100, 13 + Pi, 17 + EulerGamma], 400], N[Re[Beta[(1/3) + (Catalan I), Pi + (EulerGamma I)]], 400], N[Im[Beta[(1/3) + (Catalan I), Pi + (EulerGamma I)]], 400], N[Gamma[(1/4) + Pi, (1/3)], 400], N[Gamma[101 + Pi, (31/3)], 400], N[Gamma[(1/4) + Pi, (11/3), (29/3)], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00307_polygamma(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_gamma(10 * k, (10 * k) + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00307_polygamma",
                                       "Table[N[PolyGamma[10 k, 10 k + Pi], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00308_polygamma_fixed_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_gamma(11, ef::tenth() + (((20 * k) - 1000) * ef::pi()));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00308_polygamma_fixed_n",
                                       "Table[N[PolyGamma[11, (1/10) + (((20 k) - 1000) Pi)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00411_bessel_iv(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        static const e_float pi_over_ten         = ef::pi() / 10;"),
    std::string("        static const e_float thousand_plus_gamma = ef::thousand() + ef::euler_gamma();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(thousand_plus_gamma - (20 * k),"),
    std::string("                                                               pi_over_ten         + (20 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00411_bessel_iv",
                                       "Table[N[BesselI[(1000 - (20 k)) + EulerGamma, (20 k) + (Pi/10)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00412_bessel_jv(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        static const e_float pi_over_ten         = ef::pi() / 10;"),
    std::string("        static const e_float thousand_plus_gamma = ef::thousand() + ef::euler_gamma();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(thousand_plus_gamma - (20 * k),"),
    std::string("                                                               pi_over_ten         + (20 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00412_bessel_jv",
                                       "Table[N[BesselJ[(1000 - (20 k)) + EulerGamma, (20 k) + (Pi/10)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00413_bessel_kv(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        static const e_float pi_over_ten         = ef::pi() / 10;"),
    std::string("        static const e_float thousand_plus_gamma = ef::thousand() + ef::euler_gamma();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(thousand_plus_gamma - (20 * k),"),
    std::string("                                                               pi_over_ten         + (20 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00413_bessel_kv",
                                       "Table[N[BesselK[(1000 - (20 k)) + EulerGamma, (20 k) + (Pi/10)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00414_bessel_yv(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        static const e_float pi_over_ten         = ef::pi() / 10;"),
    std::string("        static const e_float thousand_plus_gamma = ef::thousand() + ef::euler_gamma();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(thousand_plus_gamma - (20 * k),"),
    std::string("                                                               pi_over_ten         + (20 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00414_bessel_yv",
                                       "Table[N[BesselY[(1000 - (20 k)) + EulerGamma, (20 k) + (Pi/10)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00415_bessel_jv_negative_v(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() - (10 * k), ef::pi() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00415_bessel_jv_negative_v",
                                       "Table[N[BesselJ[EulerGamma - (10 k), Pi + (10 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00416_bessel_yv_negative_v(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() - (10 * k), ef::pi() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00416_bessel_yv_negative_v",
                                       "Table[N[BesselY[EulerGamma - (10 k), Pi + (10 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00421_bessel_iv_medium_v(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + (k * 2000),"),
    std::string("                                                               ef::hundred() + (500 * (k * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00421_bessel_iv_medium_v",
                                       "Table[N[BesselI[(k 2000) + EulerGamma, 100 + (500 (k k))], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00422_bessel_jv_medium_v(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + (k * 2000),"),
    std::string("                                                               ef::hundred() + (500 * (k * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00422_bessel_jv_medium_v",
                                       "Table[N[BesselJ[(k 2000) + EulerGamma, 100 + (500 (k k))], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00423_bessel_kv_medium_v(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + (k * 2000),"),
    std::string("                                                               ef::hundred() + (500 * (k * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00423_bessel_kv_medium_v",
                                       "Table[N[BesselK[(k 2000) + EulerGamma, 100 + (500 (k k))], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00424_bessel_yv_medium_v(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + (k * 2000),"),
    std::string("                                                               ef::hundred() + (500 * (k * k)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00424_bessel_yv_medium_v",
                                       "Table[N[BesselY[(k 2000) + EulerGamma, 100 + (500 (k k))], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00431_bessel_iv_small_v_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ratl_num = (e_float(83) * k) / static_cast<INT32>(863);"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::third() + ratl_num,"),
    std::string("                                                               ef::euler_gamma() + (2 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00431_bessel_iv_small_v_small_x",
                                       "Table[N[BesselI[(1/3) + ((83 k) / 863), EulerGamma + (2 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00432_bessel_jv_small_v_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ratl_num = (e_float(83) * k) / static_cast<INT32>(863);"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::third() + ratl_num,"),
    std::string("                                                               ef::euler_gamma() + (2 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00432_bessel_jv_small_v_small_x",
                                       "Table[N[BesselJ[(1/3) + ((83 k) / 863), EulerGamma + (2 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00433_bessel_kv_small_v_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ratl_num = (e_float(83) * k) / static_cast<INT32>(863);"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::third() + ratl_num,"),
    std::string("                                                               ef::euler_gamma() + (2 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00433_bessel_kv_small_v_small_x",
                                       "Table[N[BesselK[(1/3) + ((83 k) / 863), EulerGamma + (2 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00434_bessel_yv_small_v_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ratl_num = (e_float(83) * k) / static_cast<INT32>(863);"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::third() + ratl_num,"),
    std::string("                                                               ef::euler_gamma() + (2 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00434_bessel_yv_small_v_small_x",
                                       "Table[N[BesselY[(1/3) + ((83 k) / 863), EulerGamma + (2 k)], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00441_bessel_jv_zeros(void)
{
  static const std::tr1::array<std::string, 2u> cpp_strings =
  {
    std::string("        const std::deque<e_float> zeros = ef::cyl_bessel_j_zero(ef::third() + 50, 10u);"),
    std::string("        data.assign(zeros.begin(), zeros.end());"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00441_bessel_jv_zeros",
                                       "Table[N[BesselJZero[50 + (1/3), k], 400], {k, 1, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00442_bessel_jn_zeros(void)
{
  static const std::tr1::array<std::string, 2u> cpp_strings =
  {
    std::string("        const std::deque<e_float> zeros = ef::cyl_bessel_j_zero(51, 10u);"),
    std::string("        data.assign(zeros.begin(), zeros.end());"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00442_bessel_jn_zeros",
                                       "Table[N[BesselJZero[51, k], 400], {k, 1, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00443_bessel_jv_small_v_zeros(void)
{
  static const std::tr1::array<std::string, 2u> cpp_strings =
  {
    std::string("        const std::deque<e_float> zeros = ef::cyl_bessel_j_zero(ef::tenth(), 10u);"),
    std::string("        data.assign(zeros.begin(), zeros.end());"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00443_bessel_jv_small_v_zeros",
                                       "Table[N[BesselJZero[1/10, k], 400], {k, 1, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00451_bessel_in(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        static const e_float pi_over_ten = ef::pi() / 10;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(1000 - (10 * k), pi_over_ten + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00451_bessel_in",
                                       "Table[N[BesselI[1000 - (10 k), (Pi/10) + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00452_bessel_jn(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        static const e_float pi_over_ten = ef::pi() / 10;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(1000 - (10 * k), pi_over_ten + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00452_bessel_jn",
                                       "Table[N[BesselJ[1000 - (10 k), (Pi/10) + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00453_bessel_kn(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        static const e_float pi_over_ten = ef::pi() / 10;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(500 - (5 * k), pi_over_ten + (5 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00453_bessel_kn",
                                       "Table[N[BesselK[500 - (5 k), (Pi/10) + (5 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00454_bessel_yn(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        static const e_float pi_over_ten = ef::pi() / 10;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(500 - (5 * k), pi_over_ten + (5 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00454_bessel_yn",
                                       "Table[N[BesselY[500 - (5 k), (Pi/10) + (5 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00460_bessel_jv_yv_asymp(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(10u);"),
    std::string("        INT32 n = static_cast<INT32>(1);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = n + ef::euler_gamma();"),
    std::string("          const e_float x = n + ef::pi();"),
    std::string("          const e_float v_plus_one = v + 1;"),
    std::string("          data[static_cast<std::size_t>(k)] =   (ef::cyl_bessel_j(v_plus_one, x) * ef::cyl_bessel_y(v, x))"),
    std::string("                                              - (ef::cyl_bessel_j(v, x)          * ef::cyl_bessel_y(v_plus_one, x));"),
    std::string("          n *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00460_bessel_jv_yv_asymp",
                                       "Table[N[2 / (Pi ((10^k) + Pi)), 400], {k, 0, 9, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00461_bessel_jv_asymp(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(7u);"),
    std::string("        static const e_float v = ef::ten_k() + ef::euler_gamma();"),
    std::string("        INT32 n = static_cast<INT32>(7000);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(v, n + ef::pi());"),
    std::string("          n += static_cast<INT32>(1000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00461_bessel_jv_asymp",
                                       "Table[N[BesselJ[10000 + EulerGamma, 7000 + (k * 1000) + Pi], 400], {k, 0, 6, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00462_bessel_yv_asymp(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(7u);"),
    std::string("        static const e_float v = ef::ten_k() + ef::euler_gamma();"),
    std::string("        INT32 n = static_cast<INT32>(7000);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(v, n + ef::pi());"),
    std::string("          n += static_cast<INT32>(1000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00462_bessel_yv_asymp",
                                       "Table[N[BesselY[10000 + EulerGamma, 7000 + (k * 1000) + Pi], 400], {k, 0, 6, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00463_bessel_iv_kv_asymp(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(9u);"),
    std::string("        INT32 n = static_cast<INT32>(1);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = n + ef::euler_gamma();"),
    std::string("          const e_float x = n + ef::pi();"),
    std::string("          const e_float v_plus_one = v + 1;"),
    std::string("          data[static_cast<std::size_t>(k)] =   (ef::cyl_bessel_k(v_plus_one, x) * ef::cyl_bessel_i(v, x))"),
    std::string("                                              + (ef::cyl_bessel_k(v, x)          * ef::cyl_bessel_i(v_plus_one, x));"),
    std::string("          n *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "00463_bessel_iv_kv_asymp",
                                       "Table[N[1 / ((10^k) + Pi), 400], {k, 0, 8, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00471_bessel_jn_asymp(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(7u);"),
    std::string("        static const e_float v = ef::ten_k() + ef::euler_gamma();"),
    std::string("        INT32 n = static_cast<INT32>(7000);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(static_cast<INT32>(10000), n + ef::pi());"),
    std::string("          n += static_cast<INT32>(1000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00471_bessel_jn_asymp",
                                       "Table[N[BesselJ[10000, 7000 + (k * 1000) + Pi], 400], {k, 0, 6, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00472_bessel_yn_asymp(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(7u);"),
    std::string("        static const e_float v = ef::ten_k() + ef::euler_gamma();"),
    std::string("        INT32 n = static_cast<INT32>(7000);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(static_cast<INT32>(10000), n + ef::pi());"),
    std::string("          n += static_cast<INT32>(1000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00472_bessel_yn_asymp",
                                       "Table[N[BesselY[10000, 7000 + (k * 1000) + Pi], 400], {k, 0, 6, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00481_bessel_iv_at_v_1_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00481_bessel_iv_at_v_1_3",
                                       "Table[N[BesselI[1/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00482_bessel_jv_at_v_1_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00482_bessel_jv_at_v_1_3",
                                       "Table[N[BesselJ[1/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00483_bessel_kv_at_v_1_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00483_bessel_kv_at_v_1_3",
                                       "Table[N[BesselK[1/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00484_bessel_yv_at_v_1_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00484_bessel_yv_at_v_1_3",
                                       "Table[N[BesselY[1/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00486_bessel_iv_at_v_301_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(100 + ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00486_bessel_iv_at_v_301_3",
                                       "Table[N[BesselI[301/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00487_bessel_jv_at_v_301_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(100 + ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00487_bessel_jv_at_v_301_3",
                                       "Table[N[BesselJ[301/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00488_bessel_kv_at_v_301_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(100 + ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00488_bessel_kv_at_v_301_3",
                                       "Table[N[BesselK[301/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00489_bessel_yv_at_v_301_3(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(100 + ef::third(), ef::euler_gamma() + (10 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00489_bessel_yv_at_v_301_3",
                                       "Table[N[BesselY[301/3, EulerGamma + (10 k)], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00491_bessel_iv_large_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00491_bessel_iv_large_x",
                                       "Table[N[BesselI[301/3, EulerGamma + (10^6) (k + 1)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00492_bessel_jv_large_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00492_bessel_jv_large_x",
                                       "Table[N[BesselJ[301/3, EulerGamma + (10^6) (k + 1)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00493_bessel_kv_large_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00493_bessel_kv_large_x",
                                       "Table[N[BesselK[301/3, EulerGamma + (10^6) (k + 1)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00494_bessel_yv_large_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00494_bessel_yv_large_x",
                                       "Table[N[BesselY[301/3, EulerGamma + (10^6) (k + 1)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00501_bessel_iv_order0_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + 1, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00501_bessel_iv_order0_v",
                                       "Table[N[BesselI[EulerGamma + (10^0), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00502_bessel_jv_order0_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 1, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00502_bessel_jv_order0_v",
                                       "Table[N[BesselJ[EulerGamma + (10^0), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00503_bessel_kv_order0_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + 1, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00503_bessel_kv_order0_v",
                                       "Table[N[BesselK[EulerGamma + (10^0), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00504_bessel_yv_order0_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + 1, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00504_bessel_yv_order0_v",
                                       "Table[N[BesselY[EulerGamma + (10^0), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00511_bessel_iv_order1_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + 10, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00511_bessel_iv_order1_v",
                                       "Table[N[BesselI[EulerGamma + (10^1), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00512_bessel_jv_order1_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 10, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00512_bessel_jv_order1_v",
                                       "Table[N[BesselJ[EulerGamma + (10^1), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00513_bessel_kv_order1_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + 10, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00513_bessel_kv_order1_v",
                                       "Table[N[BesselK[EulerGamma + (10^1), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00514_bessel_yv_order1_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + 10, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00514_bessel_yv_order1_v",
                                       "Table[N[BesselY[EulerGamma + (10^1), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00521_bessel_iv_order2_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + 100, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00521_bessel_iv_order2_v",
                                       "Table[N[BesselI[EulerGamma + (10^2), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00522_bessel_jv_order2_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 100, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00522_bessel_jv_order2_v",
                                       "Table[N[BesselJ[EulerGamma + (10^2), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00523_bessel_kv_order2_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + 100, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00523_bessel_kv_order2_v",
                                       "Table[N[BesselK[EulerGamma + (10^2), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00524_bessel_yv_order2_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + 100, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00524_bessel_yv_order2_v",
                                       "Table[N[BesselY[EulerGamma + (10^2), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00531_bessel_iv_order3_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + 1000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00531_bessel_iv_order3_v",
                                       "Table[N[BesselI[EulerGamma + (10^3), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00532_bessel_jv_order3_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 1000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00532_bessel_jv_order3_v",
                                       "Table[N[BesselJ[EulerGamma + (10^3), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00533_bessel_kv_order3_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + 1000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00533_bessel_kv_order3_v",
                                       "Table[N[BesselK[EulerGamma + (10^3), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00534_bessel_yv_order3_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + 1000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00534_bessel_yv_order3_v",
                                       "Table[N[BesselY[EulerGamma + (10^3), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00541_bessel_iv_order4_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_i(ef::euler_gamma() + 10000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00541_bessel_iv_order4_v",
                                       "Table[N[BesselI[EulerGamma + (10^4), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00542_bessel_jv_order4_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 10000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00542_bessel_jv_order4_v",
                                       "Table[N[BesselJ[EulerGamma + (10^4), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00543_bessel_kv_order4_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_k(ef::euler_gamma() + 10000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00543_bessel_kv_order4_v",
                                       "Table[N[BesselK[EulerGamma + (10^4), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00544_bessel_yv_order4_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k2 = k * k;"),
    std::string("          const INT32 k3 = k * k2;"),
    std::string("          const INT32 k7 = (k2 * k2) * k3;"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(ef::euler_gamma() + 10000, k7 + ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00544_bessel_yv_order4_v",
                                       "Table[N[BesselY[EulerGamma + (10^4), (k^7) + Pi], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00551_various_bessel(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::cyl_bessel_i(16600, e_float(1000) + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::cyl_bessel_i(16600, e_float(2000) + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::cyl_bessel_i(16600, e_float(4000) + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::cyl_bessel_k(ef::half() + e_float(\"1e-20\"), 12 + ef::pi()));"),
    std::string("        data.push_back(ef::cyl_bessel_k(ef::half() + e_float(\"1e-50\"), 12 + ef::pi()));"),
    std::string("        data.push_back(ef::cyl_bessel_y(123, ef::million()));"),
    std::string("        data.push_back(ef::cyl_bessel_y( 79, ef::hundred_k()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00551_various_bessel",
                                       "{N[BesselI[16600, 1000 + EulerGamma], 400], N[BesselI[16600, 2000 + EulerGamma], 400], N[BesselI[16600, 4000 + EulerGamma], 400], N[BesselK[(1/2) + (10^(-20)), 12 + Pi], 400], N[BesselK[(1/2) + (10^(-50)), 12 + Pi], 400], N[BesselY[123, 1000000], 400], N[BesselY[79, 100000], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00571_airy_ai(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_a(((20 * k) - 1000) * ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00571_airy_ai",
                                       "Table[N[AiryAi[((20 k) - 1000) Pi], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00572_airy_ai_prime(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_a_prime(((20 * k) - 1000) * ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00572_airy_ai_prime",
                                       "Table[N[AiryAiPrime[((20 k) - 1000) Pi], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00573_airy_bi(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(61u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_b(((20 * k) - 1000) * ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00573_airy_bi",
                                       "Table[N[AiryBi[((20 k) - 1000) Pi], 400], {k, 0, 60, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00574_airy_bi_prime(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(61u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_b_prime(((20 * k) - 1000) * ef::pi());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00574_airy_bi_prime",
                                       "Table[N[AiryBiPrime[((20 k) - 1000) Pi], 400], {k, 0, 60, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00581_airy_ai_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_a(((e_float(k) / 5) - 10) * ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00581_airy_ai_small_x",
                                       "Table[N[AiryAi[((k / 5) - 10) EulerGamma], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00582_airy_ai_prime_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_a_prime(((e_float(k) / 5) - 10) * ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00582_airy_ai_prime_small_x",
                                       "Table[N[AiryAiPrime[((k / 5) - 10) EulerGamma], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00583_airy_bi_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_b(((e_float(k) / 5) - 10) * ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00583_airy_bi_small_x",
                                       "Table[N[AiryBi[((k / 5) - 10) EulerGamma], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00584_airy_bi_prime_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::airy_b_prime(((e_float(k) / 5) - 10) * ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00584_airy_bi_prime_small_x",
                                       "Table[N[AiryBiPrime[((k / 5) - 10) EulerGamma], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00591_airy_ai_zeros(void)
{
  static const std::tr1::array<std::string, 2u> cpp_strings =
  {
    std::string("        const std::deque<e_float> zeros = ef::airy_a_zero(21u);"),
    std::string("        data.assign(zeros.begin(), zeros.end());"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00591_airy_ai_zeros",
                                       "Table[N[AiryAiZero[k + 1], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00592_airy_bi_zeros(void)
{
  static const std::tr1::array<std::string, 2u> cpp_strings =
  {
    std::string("        const std::deque<e_float> zeros = ef::airy_b_zero(21u);"),
    std::string("        data.assign(zeros.begin(), zeros.end());"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00592_airy_bi_zeros",
                                       "Table[N[AiryBiZero[k + 1], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00601_legendre_pvu_vary_01(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00601_legendre_pvu_vary_01",
                                       "Table[N[LegendreP[+((10 - k) * 13) + Sqrt[1/3], +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00602_legendre_pvu_vary_02(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00602_legendre_pvu_vary_02",
                                       "Table[N[LegendreP[-((10 - k) * 13) - Sqrt[1/3], -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00603_legendre_pvu_vary_03(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00603_legendre_pvu_vary_03",
                                       "Table[N[LegendreP[+((10 - k) * 13) + Sqrt[1/3], -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00604_legendre_pvu_vary_04(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00604_legendre_pvu_vary_04",
                                       "Table[N[LegendreP[-((10 - k) * 13) - Sqrt[1/3], +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00605_legendre_pvu_vary_05(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00605_legendre_pvu_vary_05",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00606_legendre_pvu_vary_06(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00606_legendre_pvu_vary_06",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00607_legendre_pvu_vary_07(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00607_legendre_pvu_vary_07",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00608_legendre_pvu_vary_08(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00608_legendre_pvu_vary_08",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00609_legendre_pvu_vary_09(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00609_legendre_pvu_vary_09",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00610_legendre_pvu_vary_10(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = -17  - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00610_legendre_pvu_vary_10",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], -17 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00611_legendre_pvu_vary_11(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = -17  - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00611_legendre_pvu_vary_11",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], -17 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00612_legendre_pvu_vary_12(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00612_legendre_pvu_vary_12",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00613_legendre_pvu_vary_13(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00613_legendre_pvu_vary_13",
                                       "Table[N[LegendreP[+17 + Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00614_legendre_pvu_vary_14(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00614_legendre_pvu_vary_14",
                                       "Table[N[LegendreP[-17 - Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00615_legendre_pvu_vary_15(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00615_legendre_pvu_vary_15",
                                       "Table[N[LegendreP[+17 + Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00616_legendre_pvu_vary_16(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00616_legendre_pvu_vary_16",
                                       "Table[N[LegendreP[-17 - Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00621_legendre_qvu_vary_01(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00621_legendre_qvu_vary_01",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13) + Sqrt[1/3], +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00622_legendre_qvu_vary_02(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00622_legendre_qvu_vary_02",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13) - Sqrt[1/3], -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00623_legendre_qvu_vary_03(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00623_legendre_qvu_vary_03",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13) + Sqrt[1/3], -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00624_legendre_qvu_vary_04(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00624_legendre_qvu_vary_04",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13) - Sqrt[1/3], +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00625_legendre_qvu_vary_05(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00625_legendre_qvu_vary_05",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00626_legendre_qvu_vary_06(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00626_legendre_qvu_vary_06",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00627_legendre_qvu_vary_07(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00627_legendre_qvu_vary_07",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00628_legendre_qvu_vary_08(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00628_legendre_qvu_vary_08",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00629_legendre_qvu_vary_09(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00629_legendre_qvu_vary_09",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00630_legendre_qvu_vary_10(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = -17  - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00630_legendre_qvu_vary_10",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], -17 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00631_legendre_qvu_vary_11(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const e_float u  = -17  - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00631_legendre_qvu_vary_11",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], -17 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00632_legendre_qvu_vary_12(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00632_legendre_qvu_vary_12",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00633_legendre_qvu_vary_13(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00633_legendre_qvu_vary_13",
                                       "Table[N[Re[LegendreQ[+17 + Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00634_legendre_qvu_vary_14(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00634_legendre_qvu_vary_14",
                                       "Table[N[Re[LegendreQ[-17 - Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00635_legendre_qvu_vary_15(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const e_float u  = -127 - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00635_legendre_qvu_vary_15",
                                       "Table[N[Re[LegendreQ[+17 + Sqrt[1/3], -127 - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00636_legendre_qvu_vary_16(void)
{
  static const std::tr1::array<std::string, 12u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(v, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00636_legendre_qvu_vary_16",
                                       "Table[N[Re[LegendreQ[-17 - Sqrt[1/3], +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00641_legendre_pvm_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00641_legendre_pvm_vary_01",
                                       "Table[N[LegendreP[+((10 - k) * 13) + Sqrt[1/3], +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00642_legendre_pvm_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const INT32   m  = -((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00642_legendre_pvm_vary_02",
                                       "Table[N[LegendreP[-((10 - k) * 13) - Sqrt[1/3], -((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00644_legendre_pvm_vary_04(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00644_legendre_pvm_vary_04",
                                       "Table[N[LegendreP[-((10 - k) * 13) - Sqrt[1/3], +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00645_legendre_pvm_vary_05(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00645_legendre_pvm_vary_05",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00648_legendre_pvm_vary_08(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00648_legendre_pvm_vary_08",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00649_legendre_pvm_vary_09(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const INT32   m  = +17;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00649_legendre_pvm_vary_09",
                                       "Table[N[LegendreP[+127 + Sqrt[1/3], +17, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00652_legendre_pvm_vary_12(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const INT32   m  = +17;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00652_legendre_pvm_vary_12",
                                       "Table[N[LegendreP[-127 - Sqrt[1/3], +17, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00653_legendre_pvm_vary_13(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00653_legendre_pvm_vary_13",
                                       "Table[N[LegendreP[+17 + Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00656_legendre_pvm_vary_16(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_p(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00656_legendre_pvm_vary_16",
                                       "Table[N[LegendreP[-17 - Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00661_legendre_qvm_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00661_legendre_qvm_vary_01",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13) + Sqrt[1/3], +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00662_legendre_qvm_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const INT32   m  = -((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00662_legendre_qvm_vary_02",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13) - Sqrt[1/3], -((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00664_legendre_qvm_vary_04(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00664_legendre_qvm_vary_04",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13) - Sqrt[1/3], +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00665_legendre_qvm_vary_05(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00665_legendre_qvm_vary_05",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00668_legendre_qvm_vary_08(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00668_legendre_qvm_vary_08",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00669_legendre_qvm_vary_09(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +127 + sqrt_1_3;"),
    std::string("          const INT32   m  = +17;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00669_legendre_qvm_vary_09",
                                       "Table[N[Re[LegendreQ[+127 + Sqrt[1/3], +17, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00672_legendre_qvm_vary_12(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -127 - sqrt_1_3;"),
    std::string("          const INT32   m  = +17;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00672_legendre_qvm_vary_12",
                                       "Table[N[Re[LegendreQ[-127 - Sqrt[1/3], +17, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00673_legendre_qvm_vary_13(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +17  + sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00673_legendre_qvm_vary_13",
                                       "Table[N[Re[LegendreQ[+17 + Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00676_legendre_qvm_vary_16(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -17  - sqrt_1_3;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_q(v, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00676_legendre_qvm_vary_16",
                                       "Table[N[Re[LegendreQ[-17 - Sqrt[1/3], +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00681_legendre_pnu_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00681_legendre_pnu_vary_01",
                                       "Table[N[LegendreP[+((10 - k) * 13), +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00682_legendre_pnu_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00682_legendre_pnu_vary_02",
                                       "Table[N[LegendreP[-((10 - k) * 13), -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00683_legendre_pnu_vary_03(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00683_legendre_pnu_vary_03",
                                       "Table[N[LegendreP[+((10 - k) * 13), -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00684_legendre_pnu_vary_04(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00684_legendre_pnu_vary_04",
                                       "Table[N[LegendreP[-((10 - k) * 13), +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00685_legendre_pnu_vary_05(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +127;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00685_legendre_pnu_vary_05",
                                       "Table[N[LegendreP[+127, +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00689_legendre_pnu_vary_09(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +127;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_p(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00689_legendre_pnu_vary_09",
                                       "Table[N[LegendreP[+127, +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00691_legendre_qnu_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00691_legendre_qnu_vary_01",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13), +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00692_legendre_qnu_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00692_legendre_qnu_vary_02",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13), -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00693_legendre_qnu_vary_03(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const e_float u  = -((10 - k) * 11) - sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00693_legendre_qnu_vary_03",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13), -((10 - k) * 11) - Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00694_legendre_qnu_vary_04(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const e_float u  = +((10 - k) * 11) + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00694_legendre_qnu_vary_04",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13), +((10 - k) * 11) + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00695_legendre_qnu_vary_05(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +127;"),
    std::string("          const e_float u  = +127 + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00695_legendre_qnu_vary_05",
                                       "Table[N[Re[LegendreQ[+127, +127 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00699_legendre_qnu_vary_09(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +127;"),
    std::string("          const e_float u  = +17  + sqrt_1_5;"),
    std::string("          data[k] = ef::legendre_q(n, u, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00699_legendre_qnu_vary_09",
                                       "Table[N[Re[LegendreQ[+127, +17 + Sqrt[1/5], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00701_legendre_pnm_vary_01(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00701_legendre_pnm_vary_01",
                                       "Table[N[LegendreP[+((10 - k) * 13), +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00702_legendre_pnm_vary_02(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const INT32   m  = -((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00702_legendre_pnm_vary_02",
                                       "Table[N[LegendreP[-((10 - k) * 13), -((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00703_legendre_pnm_vary_03(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +((10 - k) * 13);"),
    std::string("          const INT32   m  = -((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00703_legendre_pnm_vary_03",
                                       "Table[N[LegendreP[+((10 - k) * 13), -((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00704_legendre_pnm_vary_04(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = -((10 - k) * 13);"),
    std::string("          const INT32   m  = +((10 - k) * 11);"),
    std::string("          data[k] = ef::legendre_p(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00704_legendre_pnm_vary_04",
                                       "Table[N[LegendreP[-((10 - k) * 13), +((10 - k) * 11), -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00709_legendre_pnm_vary_09(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +127;"),
    std::string("          const INT32   m  = +17;"),
    std::string("          data[k] = ef::legendre_p(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00709_legendre_pnm_vary_09",
                                       "Table[N[LegendreP[+127, +17, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00721_legendre_qnm_vary_01(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x  = +(k + ef::euler_gamma()) / 11;"),
    std::string("          const INT32   n  = +((k + 1) * 4);"),
    std::string("          const INT32   m  = +((k + 1) * 3);"),
    std::string("          data[k] = ef::legendre_q(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00721_legendre_qnm_vary_01",
                                       "Table[N[LegendreQ[+((k + 1)*4), +((k + 1)*3), (k + EulerGamma)/11], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00729_legendre_qnm_vary_09(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x  = +(k + ef::euler_gamma()) / 11;"),
    std::string("          const INT32   n  = +61;"),
    std::string("          const INT32   m  = +13;"),
    std::string("          data[k] = ef::legendre_q(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00729_legendre_qnm_vary_09",
                                       "Table[N[LegendreQ[+61, +13, (k + EulerGamma)/11], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00733_legendre_qnm_vary_13(void)
{
  static const std::tr1::array<std::string, 10u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const INT32   n  = +17;"),
    std::string("          const INT32   m  = +127;"),
    std::string("          data[k] = ef::legendre_q(n, m, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00733_legendre_qnm_vary_13",
                                       "Table[N[Re[LegendreQ[+17, +127, -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00741_legendre_pv_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          data[k] = ef::legendre_p(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00741_legendre_pv_vary_01",
                                       "Table[N[LegendreP[+((10 - k) * 13) + Sqrt[1/3], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00742_legendre_pv_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          data[k] = ef::legendre_p(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00742_legendre_pv_vary_02",
                                       "Table[N[LegendreP[-((10 - k) * 13) - Sqrt[1/3], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00751_legendre_qv_vary_01(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = +((10 - k) * 13) + sqrt_1_3;"),
    std::string("          data[k] = ef::legendre_q(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00751_legendre_qv_vary_01",
                                       "Table[N[Re[LegendreQ[+((10 - k) * 13) + Sqrt[1/3], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00752_legendre_qv_vary_02(void)
{
  static const std::tr1::array<std::string, 11u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const e_float delta = ef::euler_gamma() / 100;"),
    std::string("        static const e_float sqrt_1_3 = ef::sqrt(ef::third());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float xk = ef::one_minus() +  (e_float(k) / 5);"),
    std::string("          const e_float x  = k <= 5 ? xk + delta : xk - delta;"),
    std::string("          const e_float v  = -((10 - k) * 13) - sqrt_1_3;"),
    std::string("          data[k] = ef::legendre_q(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00752_legendre_qv_vary_02",
                                       "Table[N[Re[LegendreQ[-((10 - k) * 13) - Sqrt[1/3], -1 + (k / 5) + Which[k < 6,EulerGamma/100, True, -EulerGamma/100]]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00761_various_legendre_pnm_qnm(void)
{
  static const std::tr1::array<std::string, 17u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::legendre_p(0, 0, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_p(0, 1, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_p(1, 0, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_p(1, 1, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_q(0, 0, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_q(0, 1, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_q(1, 0, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_q(1, 1, ef::euler_gamma()));"),
    std::string("        data.push_back(ef::legendre_p( 10, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p( 11, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p(500, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p(501, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_q( 11, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_q( 12, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_q(501, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_q(502, ef::zero()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00761_various_legendre_pnm_qnm",
                                       "{N[LegendreP[0, 0, EulerGamma], 400], N[LegendreP[0, 1, EulerGamma], 400], N[LegendreP[1, 0, EulerGamma], 400], N[LegendreP[1, 1, EulerGamma], 400], N[LegendreQ[0, 0, EulerGamma], 400], N[LegendreQ[0, 1, EulerGamma], 400], N[LegendreQ[1, 0, EulerGamma], 400], N[LegendreQ[1, 1, EulerGamma], 400], N[LegendreP[10, 0], 400], N[LegendreP[11, 0], 400], N[LegendreP[500, 0], 400], N[LegendreP[501, 0], 400], N[LegendreQ[11, 0], 400], N[LegendreQ[12, 0], 400], N[LegendreQ[501, 0], 400], N[LegendreQ[502, 0], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00771_legendre_pvu_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          const e_float u = ef::euler_gamma() * (k + 1);"),
    std::string("          data[k] = ef::legendre_p(v * v, u * u, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00771_legendre_pvu_zero_x",
                                       "Table[N[LegendreP[((Pi / 4) (k + 1))^2, (EulerGamma (k + 1))^2, 0], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00772_legendre_pvm_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          const INT32   m = k + 1;"),
    std::string("          data[k] = ef::legendre_p(v * v, m * m, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00772_legendre_pvm_zero_x",
                                       "Table[N[LegendreP[((Pi / 4) (k + 1))^2, (k + 1)^2, 0], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00773_legendre_pnm_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 n = static_cast<INT32>(static_cast<INT32>(2 * k) + static_cast<INT32>(1));"),
    std::string("          const INT32 m = static_cast<INT32>(k + 1);"),
    std::string("          data[k] = ef::legendre_p(n * n, m * m, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00773_legendre_pnm_zero_x",
                                       "Table[N[LegendreP[((2 k) + 1)^2, (k + 1)^2, 0], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00781_legendre_qvu_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          const e_float u = ef::euler_gamma() * (k + 1);"),
    std::string("          data[k] = ef::legendre_q(v, u, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00781_legendre_qvu_zero_x",
                                       "Table[N[Re[LegendreQ[(Pi/4) (k + 1), EulerGamma (k + 1), 0]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00782_legendre_qvm_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          const INT32   m = k + 1;"),
    std::string("          data[k] = ef::legendre_q(v, m, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00782_legendre_qvm_zero_x",
                                       "Table[N[Re[LegendreQ[(Pi / 4) (k + 1), k + 1, 0]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00783_legendre_qnm_zero_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 n = static_cast<INT32>(static_cast<INT32>(2 * k) + static_cast<INT32>(1));"),
    std::string("          const INT32 m = static_cast<INT32>(k + 1);"),
    std::string("          data[k] = ef::legendre_q(n * n, m * m, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00783_legendre_qnm_zero_x",
                                       "Table[N[Re[LegendreQ[((2 k) + 1)^2, (k + 1)^2, 0]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00788_legendre_qvm_tiny_x(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        e_float ten_pow_nk = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = ef::pi() + 31;"),
    std::string("          const INT32   m = static_cast<INT32>(39);"),
    std::string("          data[k] = ef::legendre_q(v, m, ef::half() * ten_pow_nk);"),
    std::string("          ten_pow_nk /= static_cast<INT32>(10000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00788_legendre_qvm_tiny_x",
                                       "Table[N[Re[LegendreQ[Pi + 31, 39, (1/2) / (10^(4 k))]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00791_legendre_pv_zero_x(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          data[k] = ef::legendre_p(v * v, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00791_legendre_pv_zero_x",
                                       "Table[N[LegendreP[((Pi / 4) (k + 1))^2, 0], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00792_legendre_qv_zero_x(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = (ef::pi() / 4) * (k + 1);"),
    std::string("          data[k] = ef::legendre_q(v, ef::zero());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00792_legendre_qv_zero_x",
                                       "Table[N[Re[LegendreQ[(Pi/4) (k + 1), 0]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00801_conf_hyperg_vary1_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = ef::third()   + static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = ef::quarter() + static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00801_conf_hyperg_vary1_mpx",
                                       "Table[N[Hypergeometric1F1[(1/3) + (100 - (10 k)), (1/4) + (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00802_conf_hyperg_vary2_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = -ef::quarter() - static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00802_conf_hyperg_vary2_mpx",
                                       "Table[N[Hypergeometric1F1[-(1/3) - (100 - (10 k)), -(1/4) - (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00803_conf_hyperg_vary3_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = +ef::quarter() + static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00803_conf_hyperg_vary3_mpx",
                                       "Table[N[Hypergeometric1F1[-(1/3) - (100 - (10 k)), (1/4) + (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00804_conf_hyperg_vary4_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = +ef::third()   + static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = -ef::quarter() - static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00804_conf_hyperg_vary4_mpx",
                                       "Table[N[Hypergeometric1F1[+(1/3) + (100 - (10 k)), -(1/4) - (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00811_conf_hyperg_vary1_mnx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = ef::third()   + static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = ef::quarter() + static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00811_conf_hyperg_vary1_mnx",
                                       "Table[N[Hypergeometric1F1[(1/3) + (100 - (10 k)), (1/4) + (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00821_conf_hyperg_pos_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = ef::third()   + static_cast<INT32>(50);"),
    std::string("          const e_float bk = ef::quarter() + static_cast<INT32>(50);"),
    std::string("          data[k] = ef::conf_hyperg(ak , bk, ef::euler_gamma() + (500 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00821_conf_hyperg_pos_x",
                                       "Table[N[Hypergeometric1F1[(1/3) + 50, (1/4) + 50, EulerGamma + (500 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00826_various_conf_hyperg(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::conf_hyperg(ef::zero(), ef::pi(), ef::third()));"),
    std::string("        data.push_back(ef::conf_hyperg(ef::one(), ef::pi(), ef::third()));"),
    std::string("        data.push_back(ef::conf_hyperg(ef::twenty(), ef::pi(), ef::third()));"),
    std::string("        data.push_back(ef::conf_hyperg(ef::three(), ef::pi() + (ef::five() / 4), ef::thirty() + ef::two_third()));"),
    std::string("        data.push_back(ef::conf_hyperg(ef::twenty(), ef::thirty(), ef::ten() + ef::pi()));"),
    std::string("        data.push_back(ef::conf_hyperg(ef::ten() + ef::third(), ef::ten() + ef::third(), ef::hundred() + ef::pi()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00826_various_conf_hyperg",
                                       "{N[Re[Hypergeometric1F1[0, Pi, (1/3)]], 400], N[Re[Hypergeometric1F1[1, Pi, (1/3)]], 400], N[Re[Hypergeometric1F1[20, Pi, (1/3)]], 400], N[Re[Hypergeometric1F1[3, Pi + (5/4), 30 + (2/3)]], 400], N[Re[Hypergeometric1F1[20, 30, 10 + Pi]], 400], N[Re[Hypergeometric1F1[10 + (1/3), 10 + (1/3), 100 + Pi]], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00831_laguerre_vary1_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = ef::third()   + static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = ef::quarter() + static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00831_laguerre_vary1_mpx",
                                       "Table[N[LaguerreL[(1/3) + (100 - (10 k)), (1/4) + (0 + (10 k)), EulerGamma + (100 k)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00832_laguerre_vary2_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = -ef::quarter() - static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00832_laguerre_vary2_mpx",
                                       "Table[N[Re[LaguerreL[-(1/3) - (100 - (10 k)), -(1/4) - (0 + (10 k)), EulerGamma + (100 k)]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00833_laguerre_vary3_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = +ef::quarter() + static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00833_laguerre_vary3_mpx",
                                       "Table[N[Re[LaguerreL[-(1/3) - (100 - (10 k)), (1/4) + (0 + (10 k)), EulerGamma + (100 k)]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00834_laguerre_vary4_mpx(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float ak = +ef::third()   + static_cast<INT32>(100 - (10 * k));"),
    std::string("          const e_float bk = -ef::quarter() - static_cast<INT32>(0   + (10 * k));"),
    std::string("          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00834_laguerre_vary4_mpx",
                                       "Table[N[Re[LaguerreL[+(1/3) + (100 - (10 k)), -(1/4) - (0 + (10 k)), EulerGamma + (100 k)]], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00835_various_laguerre(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::laguerre(ef::euler_gamma() + 77, ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(ef::fifty(), ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(50, 60, ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(11, 123, ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(123, 11, ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(50, 60 + ef::euler_gamma(), ef::hundred() + ef::pi()));"),
    std::string("        data.push_back(ef::laguerre(13, -23 - ef::euler_gamma(), ef::hundred() + ef::pi()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00835_various_laguerre",
                                       "{N[Re[LaguerreL[EulerGamma + 77, 100 + Pi]], 400], N[Re[LaguerreL[50, 100 + Pi]], 400], N[Re[LaguerreL[50, 60, 100 + Pi]], 400], N[Re[LaguerreL[11, 123, 100 + Pi]], 400], N[Re[LaguerreL[123, 11, 100 + Pi]], 400], N[Re[LaguerreL[50, 60 + EulerGamma, 100 + Pi]], 400], N[Re[LaguerreL[13, -23 - EulerGamma, 100 + Pi]], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00851_weber_d_pos_x_pos_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = +ef::euler_gamma() + (k * static_cast<INT32>(10));"),
    std::string("          const e_float x = +pi_over_ten       + (k * static_cast<INT32>(20));"),
    std::string("          data[k] = ef::weber_d(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00851_weber_d_pos_x_pos_v",
                                       "Table[N[ParabolicCylinderD[+EulerGamma + (k 10), +(Pi/10) + (k 20)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00852_weber_d_neg_x_neg_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = -ef::euler_gamma() - (k * static_cast<INT32>(10));"),
    std::string("          const e_float x = -pi_over_ten       - (k * static_cast<INT32>(20));"),
    std::string("          data[k] = ef::weber_d(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00852_weber_d_neg_x_neg_v",
                                       "Table[N[ParabolicCylinderD[-EulerGamma - (k 10), -(Pi/10) - (k 20)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00853_weber_d_neg_x_pos_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = +ef::euler_gamma() + (k * static_cast<INT32>(10));"),
    std::string("          const e_float x = -pi_over_ten       - (k * static_cast<INT32>(20));"),
    std::string("          data[k] = ef::weber_d(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00853_weber_d_neg_x_pos_v",
                                       "Table[N[ParabolicCylinderD[+EulerGamma + (k 10), -(Pi/10) - (k 20)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00854_weber_d_pos_x_neg_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = -ef::euler_gamma() - (k * static_cast<INT32>(10));"),
    std::string("          const e_float x = +pi_over_ten       + (k * static_cast<INT32>(20));"),
    std::string("          data[k] = ef::weber_d(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00854_weber_d_pos_x_neg_v",
                                       "Table[N[ParabolicCylinderD[-EulerGamma - (k 10), +(Pi/10) + (k 20)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00855_weber_d_med_x_med_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = static_cast<INT32>((k + 1) * 100) + ef::pi();"),
    std::string("          const e_float x = 150 + ef::third();"),
    std::string("          data[k] = ef::weber_d(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00855_weber_d_med_x_med_v",
                                       "Table[N[ParabolicCylinderD[((k + 1) 100) + Pi, 150 + (1/3)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00861_weber_d_asymp_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(6u);"),
    std::string("        INT32 p10 = static_cast<INT32>(1);"),
    std::string("        const e_float x = ef::euler_gamma() + static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[k] = ef::weber_d(+ef::pi() + p10, x);"),
    std::string("          p10 *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00861_weber_d_asymp_v",
                                       "Table[N[ParabolicCylinderD[+Pi + (10^k), EulerGamma + 10], 400], {k, 0, 5, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00871_hermite_pos_x_pos_v(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float v = +ef::euler_gamma() + (k * static_cast<INT32>(10));"),
    std::string("          const e_float x = +pi_over_ten       + (k * static_cast<INT32>(20));"),
    std::string("          data[k] = ef::hermite(v, x);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00871_hermite_pos_x_pos_v",
                                       "Table[N[HermiteH[+EulerGamma + (k 10), +(Pi/10) + (k 20)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00901_zeta_small_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::riemann_zeta((ef::pi() * (k + 1)) / 100);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00901_zeta_small_x",
                                       "Table[N[Zeta[(Pi (k + 1)) / 100], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00902_zeta_all_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::riemann_zeta((ef::pi() * k) / 2);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00902_zeta_all_x",
                                       "Table[N[Zeta[(Pi k) / 2], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00903_zeta_neg_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::riemann_zeta((-103 * e_float(k)) / 227);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00903_zeta_neg_x",
                                       "Table[N[Zeta[(-103 k) / 227], 400], {k, 0, 50, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00931_hurw_zeta_all_x_small_a(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::hurwitz_zeta(ef::pi() * ((10 * k) + 1), ef::two_third());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00931_hurw_zeta_all_x_small_a",
                                       "Table[N[Zeta[Pi ((10 k) + 1), (2/3)], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00932_hurw_zeta_all_x_large_a(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::hurwitz_zeta(ef::pi() * ((10 * k) + 1), ef::two_third() + 1000);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00932_hurw_zeta_all_x_large_a",
                                       "Table[N[Zeta[Pi ((10 k) + 1), (2/3) + 1000], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00936_various_hurw_zeta(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::hurwitz_zeta(-101 - ef::catalan(), ef::two_third() + 80));"),
    std::string("        data.push_back(ef::hurwitz_zeta(ef::catalan(), 0));"),
    std::string("        data.push_back(ef::hurwitz_zeta(ef::catalan(), 1));"),
    std::string("        data.push_back(ef::hurwitz_zeta(ef::catalan(), 2));"),
    std::string("        data.push_back(ef::hurwitz_zeta(ef::catalan(), 17));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00936_various_hurw_zeta",
                                       "{N[Zeta[-101 - Catalan, ((2/3) + 80)], 400], N[Zeta[Catalan, 0], 400], N[Zeta[Catalan, 1], 400], N[Zeta[Catalan, 2], 400], N[Zeta[Catalan, 17], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00951_polylog_med_x_med_pos_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(20, (-20 + k) + ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00951_polylog_med_x_med_pos_n",
                                       "Table[N[PolyLog[20, (-20 + k) + EulerGamma], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00952_polylog_med_x_small_pos_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(k, (ef::euler_gamma() + k) / 13);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00952_polylog_med_x_small_pos_n",
                                       "Table[N[PolyLog[k, (EulerGamma + k) / 13], 400], {k, 0, 10, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00953_polylog_neg_x_pos_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm((k * 10) + 1, -(20 * k) - ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00953_polylog_neg_x_pos_n",
                                       "Table[N[PolyLog[(k * 10) + 1, -(20 * k) - EulerGamma], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00954_polylog_neg_x_neg_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-(k * 5) - 1, -(20 * k) - ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00954_polylog_neg_x_neg_n",
                                       "Table[N[PolyLog[-(k * 5) - 1, -(20 * k) - EulerGamma], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00955_polylog_pos_x_neg_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-(k * 5) - 1, +(20 * k) + ef::euler_gamma());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00955_polylog_pos_x_neg_n",
                                       "Table[N[PolyLog[-(k * 5) - 1, +(20 * k) + EulerGamma], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00956_polylog_x_small_neg_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-((static_cast<INT32>(2) * k) + static_cast<INT32>(1)), (k + ef::euler_gamma()) / static_cast<INT32>(20));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00956_polylog_x_small_neg_n",
                                       "Table[N[PolyLog[-((2 k) + 1), (k + EulerGamma) / 20], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00957_polylog_x_large_neg_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(6u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-500 -(k * 100), ef::five_k());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00957_polylog_x_large_neg_n",
                                       "Table[N[PolyLog[-500 - (k 100), 5000], 400], {k, 0, 5, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_00958_polylog_x_one_neg_n(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(21u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-50 - static_cast<INT32>(k * 3), ef::one() + (ef::one() / static_cast<INT32>(10000)));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_00958_polylog_x_one_neg_n",
                                       "Table[N[PolyLog[-50 - (3 k), 1 + (1 / 10000)], 400], {k, 0, 20, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01001_comp_ellint_1(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(10u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::comp_ellint_1((ef::euler_gamma() + k) / 10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01001_comp_ellint_1",
                                       "Table[N[EllipticK[(EulerGamma + k) / 10], 400], {k, 0, 9, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01002_ellint_1(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(10u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::ellint_1((ef::euler_gamma() + k) / 10, ef::two_third());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01002_ellint_1",
                                       "Table[N[EllipticF[2/3, (EulerGamma + k) / 10], 400], {k, 0, 9, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01003_comp_ellint_2(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(10u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::comp_ellint_2((ef::euler_gamma() + k) / 10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01003_comp_ellint_2",
                                       "Table[N[EllipticE[(EulerGamma + k) / 10], 400], {k, 0, 9, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01004_ellint_2(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(10u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::ellint_2((ef::euler_gamma() + k) / 10, ef::two_third());"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01004_ellint_2",
                                       "Table[N[EllipticE[2/3, (EulerGamma + k) / 10], 400], {k, 0, 9, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01021_poly_legendre_p(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::legendre_p(k, (ef::euler_gamma() + k) / 103);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01021_poly_legendre_p",
                                       "Table[N[LegendreP[k, (EulerGamma + k) / 103], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01022_poly_legendre_q(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::legendre_q(k, (ef::euler_gamma() + k) / 103);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01022_poly_legendre_q",
                                       "Table[N[LegendreQ[k, (EulerGamma + k) / 103], 400], {k, 0, 100, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01023_poly_chebyshev_t(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::chebyshev_t(k, ef::euler_gamma() + (k - 10));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01023_poly_chebyshev_t",
                                       "Table[N[ChebyshevT[k, EulerGamma + (k - 10)], 400], {k, 0, 110, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01024_poly_chebyshev_u(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::chebyshev_u(k, ef::euler_gamma() + (k - 10));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01024_poly_chebyshev_u",
                                       "Table[N[ChebyshevU[k, EulerGamma + (k - 10)], 400], {k, 0, 110, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01025_poly_hermite_h(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::hermite(k, ef::euler_gamma() + (k - 10));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01025_poly_hermite_h",
                                       "Table[N[HermiteH[k, EulerGamma + (k - 10)], 400], {k, 0, 110, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01026_poly_laguerre_l(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = ef::laguerre(k, ef::euler_gamma() + (k - 10));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01026_poly_laguerre_l",
                                       "Table[N[LaguerreL[k, EulerGamma + (k - 10)], 400], {k, 0, 110, 1}]");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01027_various_poly_legendre_p(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::legendre_p(30000, ef::third()));"),
    std::string("        data.push_back(ef::legendre_p(1202, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p(1203, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p(1204, ef::zero()));"),
    std::string("        data.push_back(ef::legendre_p(59, ef::million()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01027_various_poly_legendre_p",
                                       "{N[LegendreP[30000, (1/3)], 400], N[LegendreP[1202, 0], 400], N[LegendreP[1203, 0], 400], N[LegendreP[1204, 0], 400], N[LegendreP[59, 1000000], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_01028_various_poly(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef::chebyshev_t(-10, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::chebyshev_t(-11, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::chebyshev_u(-10, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::chebyshev_u(-11, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::hermite    (-10, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::hermite    (-11, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::laguerre   (-10, ef::ten() + ef::euler_gamma()));"),
    std::string("        data.push_back(ef::laguerre   (-11, ef::ten() + ef::euler_gamma()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_01028_various_poly",
                                       "{N[ChebyshevT[-10, 10 + EulerGamma], 400], N[ChebyshevT[-11, 10 + EulerGamma], 400], N[ChebyshevU[-10, 10 + EulerGamma], 400], N[ChebyshevU[-11, 10 + EulerGamma], 400], N[HermiteH[-10, 10 + EulerGamma], 400], N[HermiteH[-11, 10 + EulerGamma], 400], N[LaguerreL[-10, 10 + EulerGamma], 400], N[LaguerreL[-11, 10 + EulerGamma], 400]}");

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02101_z_sin(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::half() + (ef::euler_gamma() * ((100 * k) - 5000));"),
    std::string("          const e_float y = (ef::one() + (31 * k)) / static_cast<INT32>(3);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::sin(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02101_z_sin",
                                       "Table[N[Sin[((1/2) + (EulerGamma ((100 k) - 5000))) + (((1 + (31 k))/3) I)], 400], {k, 0, 100, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02102_z_cos(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::half() + (ef::euler_gamma() * ((100 * k) - 5000));"),
    std::string("          const e_float y = (ef::one() + (31 * k)) / static_cast<INT32>(3);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::cos(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02102_z_cos",
                                       "Table[N[Cos[((1/2) + (EulerGamma ((100 k) - 5000))) + (((1 + (31 k))/3) I)], 400], {k, 0, 100, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02103_z_exp(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k_plus_one = static_cast<INT32>(k + static_cast<INT32>(1));"),
    std::string("          const e_float x = ef::sqrt((ef::pi() * (100 * k_plus_one)) * (100 * k_plus_one));"),
    std::string("          const e_float y = (ef::one() + (31 * k)) / static_cast<INT32>(3);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::exp(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02103_z_exp",
                                       "Table[N[Exp[Sqrt[Pi (100 (k + 1)) (100 (k + 1))] + (((1 + (31 k)) / 3) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02104_z_log(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::tenth() + ((ef::pi() * (100 * k)) * (100 * k));"),
    std::string("          const e_float y = (ef::one() + (31 * k)) / static_cast<INT32>(3);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::log(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02104_z_log",
                                       "Table[N[Log[((1/10) + (Pi (100 k) (100 k))) + (((1 + (31 k)) / 3) I)], 400], {k, 0, 100, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02105_z_sqrt(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k_plus_one = static_cast<INT32>(k + static_cast<INT32>(1));"),
    std::string("          const e_float x = ef::pi() * k_plus_one;"),
    std::string("          const e_float y = ef::euler_gamma() * k_plus_one;"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::sqrt(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02105_z_sqrt",
                                       "Table[N[Sqrt[(Pi (k + 1)) + ((EulerGamma (k + 1)) I)], 400], {k, 0, 100, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02106_z_rootn(void)
{
  static const std::tr1::array<std::string, 8u> cpp_strings =
  {
    std::string("        data.resize(101u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const INT32 k_plus_one = static_cast<INT32>(k + 1);"),
    std::string("          const e_float x = ef::googol() * (k_plus_one * ef::pi());"),
    std::string("          const e_float y = ef::googol() * (k_plus_one * ef::euler_gamma());"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::rootn(ef_complex(x, y), k_plus_one);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02106_z_rootn",
                                       "Table[N[Power[((10^100) (Pi (k + 1))) + (((10^100) (EulerGamma (k + 1))) I), 1/(k + 1)], 400], {k, 0, 100, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02111_z_asin(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::catalan() + static_cast<INT32>(k * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::asin(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02111_z_asin",
                                       "Table[N[ArcSin[(EulerGamma + k) + ((Catalan + (k k)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02112_z_acos(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::catalan() + static_cast<INT32>(k * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::acos(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02112_z_acos",
                                       "Table[N[ArcCos[(EulerGamma + k) + ((Catalan + (k k)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02113_z_atan(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::catalan() + static_cast<INT32>(k * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::atan(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02113_z_atan",
                                       "Table[N[ArcTan[(EulerGamma + k) + ((Catalan + (k k)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02114_z_various_trig(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(efz::tan(ef_complex(ef::euler_gamma() * ef::three(), ef::golden_ratio() * ef::seven())));"),
    std::string("        data.push_back(efz::cot(ef_complex(ef::euler_gamma() * ef::three(), ef::golden_ratio() * ef::seven())));"),
    std::string("        data.push_back(efz::sec(ef_complex(ef::euler_gamma() * ef::three(), ef::golden_ratio() * ef::seven())));"),
    std::string("        data.push_back(efz::csc(ef_complex(ef::euler_gamma() * ef::three(), ef::golden_ratio() * ef::seven())));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02114_z_various_trig",
                                       "{N[Tan[(EulerGamma 3) + ((GoldenRatio 7) I)], 400], N[Cot[(EulerGamma 3) + ((GoldenRatio 7) I)], 400], N[Sec[(EulerGamma 3) + ((GoldenRatio 7) I)], 400], N[Csc[(EulerGamma 3) + ((GoldenRatio 7) I)], 400]}",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02115_z_various_elem_trans_log(void)
{
  static const std::tr1::array<std::string, 3u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(efz::log10(ef_complex(ef::euler_gamma() * ef::thousand(), ef::golden_ratio() * ef::hundred())));"),
    std::string("        data.push_back(efz::loga (ef_complex(ef::third(), ef::catalan()), ef_complex(ef::euler_gamma() * ef::thousand(), ef::golden_ratio() * ef::hundred())));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02115_z_various_elem_trans_log",
                                       "{N[Log[10, (EulerGamma 1000) + ((GoldenRatio 100) I)],400], N[Log[(1/3) + (Catalan I), (EulerGamma 1000) + ((GoldenRatio 100) I)],400]}",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02116_z_various_elem(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(ef_complex(ef::euler_gamma(), ef::golden_ratio()) + static_cast<INT32>(123));"),
    std::string("        data.push_back(ef_complex(ef::euler_gamma(), ef::golden_ratio()) - static_cast<INT32>(123));"),
    std::string("        data.push_back(static_cast<INT32>(123) + ef_complex(ef::euler_gamma(), ef::golden_ratio()));"),
    std::string("        data.push_back(static_cast<INT32>(123) - ef_complex(ef::euler_gamma(), ef::golden_ratio()));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02116_z_various_elem",
                                       "{N[(EulerGamma + (GoldenRatio I)) + 123, 400], N[(EulerGamma + (GoldenRatio I)) - 123, 400], N[123 + (EulerGamma + (GoldenRatio I)), 400], N[123 - (EulerGamma + (GoldenRatio I)), 400]}",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02121_z_sinh(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::catalan() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::sinh(ef_complex(x * x, y * y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02121_z_sinh",
                                       "Table[N[Sinh[((EulerGamma + k)^2) + (((Catalan + (2 k))^2) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02122_z_cosh(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::catalan() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::cosh(ef_complex(x * x, y * y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02122_z_cosh",
                                       "Table[N[Cosh[((EulerGamma + k)^2) + (((Catalan + (2 k))^2) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02123_z_tanh(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = (ef::euler_gamma() + k) / static_cast<INT32>(10);"),
    std::string("          const e_float y = (ef::catalan() + static_cast<INT32>(2 * k)) / static_cast<INT32>(10);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::tanh(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02123_z_tanh",
                                       "Table[N[Tanh[((EulerGamma + k) / 10) + (((Catalan + (2 k)) / 10) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02124_z_asinh(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_k = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() * ten_pow_k;"),
    std::string("          const e_float y = ef::catalan() * (ten_pow_k * static_cast<INT32>(3));"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::asinh(ef_complex(x, y));"),
    std::string("          ten_pow_k *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02124_z_asinh",
                                       "Table[N[ArcSinh[(EulerGamma 10^k) + ((Catalan 3 (10^k)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02125_z_acosh(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_k = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() * ten_pow_k;"),
    std::string("          const e_float y = ef::catalan() * (ten_pow_k * static_cast<INT32>(3));"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::acosh(ef_complex(x, y));"),
    std::string("          ten_pow_k *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02125_z_acosh",
                                       "Table[N[ArcCosh[(EulerGamma 10^k) + ((Catalan 3 (10^k)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02126_z_atanh(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        e_float ten_pow_k = ef::one();"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = (ef::euler_gamma() + k) / static_cast<INT32>(53);"),
    std::string("          const e_float y = (ef::catalan()     + k) / static_cast<INT32>(53);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::atanh(ef_complex(x, y));"),
    std::string("          ten_pow_k *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02126_z_atanh",
                                       "Table[N[ArcTanh[((EulerGamma + k) / 53) + (((Catalan + k) / 53) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02201_z_gamma(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::tenth()   + (((10 * k) - 250) * ef::pi());"),
    std::string("          const e_float y = ef::quarter() + (((10 * k) - 250) * ef::euler_gamma());"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::gamma(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02201_z_gamma",
                                       "Table[N[Gamma[((1/10) + (((10 k) - 250) Pi)) + (((1/4) + (((10 k) - 250) EulerGamma)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02202_z_gamma_medium_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(26u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::third() + (static_cast<INT32>(2) * k);"),
    std::string("          const e_float y = (ef::one() / static_cast<INT32>(7)) + (static_cast<INT32>(2) * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::gamma(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02202_z_gamma_medium_x",
                                       "Table[N[Gamma[((1/3) + (2 k)) + (((1/7) + (2 k)) I)], 400], {k, 0, 25, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02901_z_zeta_small_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = (ef::pi() * (k + 1)) / 100;"),
    std::string("          const e_float y = (ef::euler_gamma() * (k + 1)) / 100;"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::riemann_zeta(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02901_z_zeta_small_x",
                                       "Table[N[Zeta[((Pi (k + 1)) / 100) + (((EulerGamma (k + 1)) / 100) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02902_z_zeta_all_x(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(51u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::pi() * static_cast<INT32>(k + static_cast<INT32>(1));"),
    std::string("          const e_float y = ef::euler_gamma() * static_cast<INT32>(k + static_cast<INT32>(1));"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::riemann_zeta(ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02902_z_zeta_all_x",
                                       "Table[N[Zeta[(Pi (k + 1)) + ((EulerGamma (k + 1)) I)], 400], {k, 0, 50, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02903_z_zeta_neg_x(void)
{
  static const std::tr1::array<std::string, 5u> cpp_strings =
  {
    std::string("        data.resize(26u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::riemann_zeta(ef_complex((-103 * e_float(2 * k)) / 227, k + ef::euler_gamma()));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02903_z_zeta_neg_x",
                                       "Table[N[Zeta[((-103 (2 k)) / 227) + ((k + EulerGamma) I)], 400], {k, 0, 25, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02911_z_zeta_crit_strip(void)
{
  static const std::tr1::array<std::string, 9u> cpp_strings =
  {
    std::string("        data.resize(5u);"),
    std::string("        INT32 ten_pow_k = static_cast<INT32>(1);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::half();"),
    std::string("          const e_float y = ef::euler_gamma() + ten_pow_k;"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::riemann_zeta(ef_complex(x, y));"),
    std::string("          ten_pow_k *= static_cast<INT32>(10);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02911_z_zeta_crit_strip",
                                       "Table[N[Zeta[(1/2) + ((EulerGamma + 10^k) I)], 400], {k, 0, 4, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02931_z_hurw_zeta_all_x_small_a(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const ef_complex a(ef::two_third(), ef::catalan());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const ef_complex s(ef::pi() * ((10 * k) + 1), ef::euler_gamma() + k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::hurwitz_zeta(s, a);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02931_z_hurw_zeta_all_x_small_a",
                                       "Table[N[Zeta[(Pi ((10 k) + 1)) + ((EulerGamma + k) I), (2/3) + (Catalan I)], 400], {k, 0, 10, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02932_z_hurw_zeta_all_x_large_a(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(11u);"),
    std::string("        static const ef_complex a(ef::two_third() + 1000, ef::catalan());"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const ef_complex s(ef::pi() * ((10 * k) + 1), ef::euler_gamma() + k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::hurwitz_zeta(s, a);"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02932_z_hurw_zeta_all_x_large_a",
                                       "Table[N[Zeta[(Pi ((10 k) + 1)) + ((EulerGamma + k) I), ((2/3) + 1000) + (Catalan I)], 400], {k, 0, 10, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_02936_z_various_hurw_zeta(void)
{
  static const std::tr1::array<std::string, 6u> cpp_strings =
  {
    std::string("        data.clear();"),
    std::string("        data.push_back(efz::hurwitz_zeta(ef_complex(e_float(-101), -(ef::euler_gamma() + 200)), ef_complex(ef::two_third() + 80)));"),
    std::string("        data.push_back(efz::hurwitz_zeta(ef_complex(ef::catalan(), ef::euler_gamma()),  0));"),
    std::string("        data.push_back(efz::hurwitz_zeta(ef_complex(ef::catalan(), ef::euler_gamma()),  1));"),
    std::string("        data.push_back(efz::hurwitz_zeta(ef_complex(ef::catalan(), ef::euler_gamma()),  2));"),
    std::string("        data.push_back(efz::hurwitz_zeta(ef_complex(ef::catalan(), ef::euler_gamma()), 17));"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_02936_z_various_hurw_zeta",
                                       "{N[Zeta[-101 - ((EulerGamma + 200) I), ((2/3) + 80)], 400], N[Zeta[Catalan + (EulerGamma I), 0], 400], N[Zeta[Catalan + (EulerGamma I), 1], 400], N[Zeta[Catalan + (EulerGamma I), 2], 400], N[Zeta[Catalan + (EulerGamma I), 17], 400]}",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_03023_z_poly_chebyshev_t(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + (k - 10);"),
    std::string("          const e_float y = ef::golden_ratio() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::chebyshev_t(k + 1, ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_03023_z_poly_chebyshev_t",
                                       "Table[N[ChebyshevT[k + 1, (EulerGamma + (k - 10)) + ((GoldenRatio + (2 k)) I)], 400], {k, 0, 110, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_03024_z_poly_chebyshev_u(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + (k - 10);"),
    std::string("          const e_float y = ef::golden_ratio() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::chebyshev_u(k + 1, ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_03024_z_poly_chebyshev_u",
                                       "Table[N[ChebyshevU[k + 1, (EulerGamma + (k - 10)) + ((GoldenRatio + (2 k)) I)], 400], {k, 0, 110, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_03025_z_poly_hermite_h(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(111u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + (k - 10);"),
    std::string("          const e_float y = ef::golden_ratio() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::hermite(k + 1, ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_03025_z_poly_hermite_h",
                                       "Table[N[HermiteH[k + 1, (EulerGamma + (k - 10)) + ((GoldenRatio + (2 k)) I)], 400], {k, 0, 110, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

const TestCaseDescription& TestCaseDescription_case_03026_z_poly_laguerre_l(void)
{
  static const std::tr1::array<std::string, 7u> cpp_strings =
  {
    std::string("        data.resize(31u);"),
    std::string("        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)"),
    std::string("        {"),
    std::string("          const e_float x = ef::euler_gamma() + k;"),
    std::string("          const e_float y = ef::golden_ratio() + static_cast<INT32>(2 * k);"),
    std::string("          data[static_cast<std::size_t>(k)] = efz::laguerre(k + 1, ef_complex(x, y));"),
    std::string("        }"),
  };

  static const TestCaseDescription tcd(cpp_strings.begin(),
                                       cpp_strings.end(),
                                       "case_03026_z_poly_laguerre_l",
                                       "Table[N[LaguerreL[k + 1, (EulerGamma + k) + ((GoldenRatio + (2 k)) I)], 400], {k, 0, 30, 1}]",
                                       TestCaseDescription::is_complex);

  return tcd;
}

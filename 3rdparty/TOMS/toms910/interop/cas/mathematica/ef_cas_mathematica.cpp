#include <sstream>

#include <interop/cas/mathematica/ef_cas_mathematica.h>
#include <utility/util_lexical_cast.h>

void ef_cas::Mathematica::Replace_StarHat_With_e(std::string& str)
{
  static const std::string str_star_hat = "*^";
  static const std::string str_e        = "e";

  const std::string::size_type pos = str.find(str_star_hat);

  if(pos != std::string::npos)
  {
    str.replace(pos, str_star_hat.length(), str_e);
  }
}

bool ef_cas::Mathematica::SendCommand(const std::string& str_cmd) const
{
  // Send a command sequence which will command Mathematica to generate
  // a list of high precision function values.
  static_cast<void>(PutFunction("EvaluatePacket", 1));
  static_cast<void>(PutFunction("ToExpression", 1));
  static_cast<void>(PutString(str_cmd.c_str()));
  static_cast<void>(EndPacket());

  // Skip any packets before the first ReturnPacket.
  for(;;)
  {
    const int pkt = NextPacket();

    if(!pkt || pkt == ReturnPacketId())
    {
      break;
    }

    static_cast<void>(NewPacket());

    if(Error() != 0)
    {
      return false;
    }
  }

  return true;
}

void ef_cas::Mathematica::create_mp_string(std::string& str, const e_float& x) const
{
  std::stringstream ss;

  static_cast<void>(ss.precision(my_digits()));

  ss << x;

  // Create a string such as: (1.23456e49`130)
  str = "(" + ss.str() + "`" + Util::lexical_cast(my_digits()) + ")";
}

void ef_cas::Mathematica::create_mp_string(std::string& str, const ef_complex& z) const
{
  std::string str_real;
  std::string str_imag;

  create_mp_string(str_real, z.real());
  create_mp_string(str_imag, z.imag());

  // Create a string such as: ((1.23456e49`130) + ((4.5678e3`130) I))
  str = "(" + str_real + " + (" + str_imag + " I))";
}

bool ef_cas::Mathematica::get_values(const std::string& str_cmd, std::vector<e_float>& values) const
{
  values.clear();

  if(!SendCommand(str_cmd))
  {
    return false;
  }

  long len;
  if(!CheckFunction("List", len))
  {
    return false;
  }

  // Get the result strings and convert these to e_float values.
  for(std::size_t i = static_cast<std::size_t>(0u); i < static_cast<std::size_t>(len); i++)
  {
    std::string str;

    if(GetString(str))
    {
      Replace_StarHat_With_e(str);

      values.push_back(e_float(str));
    }
    else
    {
      return false;
    }
  }

  return !values.empty();
}

bool ef_cas::Mathematica::get_values(const std::string& str_cmd, std::vector<ef_complex>& values) const
{
  values.clear();

  if(!SendCommand(str_cmd))
  {
    return false;
  }

  long len;
  if(!CheckFunction("List", len))
  {
    return false;
  }

  // Get the result strings and convert these to ef_complex values.
  for(std::size_t i = static_cast<std::size_t>(0u); i < static_cast<std::size_t>(len); i++)
  {
    for(std::size_t j = static_cast<std::size_t>(0u); j < static_cast<std::size_t>(10u); j++)
    {
      const int ty = GetType();

      if(ty == 42)
      {
        break;
      }

      static_cast<void>(GetNext());
    }

    std::string str_real;
    std::string str_imag;

    if(GetString(str_real) && GetString(str_imag))
    {
      Replace_StarHat_With_e(str_real);
      Replace_StarHat_With_e(str_imag);

      values.push_back(ef_complex(e_float(str_real), e_float(str_imag)));
    }
    else
    {
      return false;
    }
  }

  return !values.empty();
}

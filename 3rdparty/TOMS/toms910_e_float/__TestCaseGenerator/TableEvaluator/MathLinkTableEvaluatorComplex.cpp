#include <MathLink/MathLinkMathematica.h>
#include <TableEvaluator/MathLinkTableEvaluatorComplex.h>

bool MathLinkTableEvaluatorComplex::Evaluate(void) const
{
  // Clear any possible values which might be present. This is necessary
  // if this expression is used for multiple data exchanges with MathLink
  // during the object lifetime.
  str_values.clear();

  // Send a command sequence which will command Mathematica to generate
  // a list of high precision function values.

  static_cast<void>(my_ML.PutFunction("EvaluatePacket", 1));
  static_cast<void>(my_ML.PutFunction("ToExpression", 1));
  static_cast<void>(my_ML.PutString(my_expression.c_str()));
  static_cast<void>(my_ML.EndPacket());

  // Skip any packets before the first ReturnPacket.
  for(;;)
  {
    const int pkt = my_ML.NextPacket();

    if(!pkt || pkt == my_ML.ReturnPacketId())
    {
      break;
    }

    static_cast<void>(my_ML.NewPacket());

    if(my_ML.Error() != 0)
    {
      return false;
    }
  }

  long len;
  if(!my_ML.CheckFunction("List", len))
  {
    return false;
  }

  str_values.resize(static_cast<std::size_t>(len));

  // Get the result strings.
  for(std::size_t i = static_cast<std::size_t>(0u); i < str_values.size(); i++)
  {
    for(std::size_t j = static_cast<std::size_t>(0u); j < static_cast<std::size_t>(10u); j++)
    {
      const int ty = my_ML.GetType();

      if(ty == 42)
      {
        break;
      }

      const int nx = my_ML.GetNext();
    }

    std::string str_real;
    std::string str_imag;

    if(my_ML.GetString(str_real) && my_ML.GetString(str_imag))
    {
      Replace_StarHat_With_e(str_real);
      Replace_StarHat_With_e(str_imag);

      str_values[i].push_back(str_real);
      str_values[i].push_back(str_imag);
    }
    else
    {
      return false;
    }
  }

  return true;
}

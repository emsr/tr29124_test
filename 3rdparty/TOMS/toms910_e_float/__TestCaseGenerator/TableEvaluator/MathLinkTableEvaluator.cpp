#include <MathLink/MathLinkMathematica.h>
#include <TableEvaluator/MathLinkTableEvaluator.h>

bool MathLinkTableEvaluator::Evaluate(void) const
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
    std::string str;

    if(my_ML.GetString(str))
    {
      Replace_StarHat_With_e(str);

      str_values[i].push_back(str);
    }
    else
    {
      return false;
    }
  }

  return true;
}

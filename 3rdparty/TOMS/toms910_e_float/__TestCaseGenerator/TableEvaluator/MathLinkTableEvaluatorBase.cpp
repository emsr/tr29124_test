#include <TableEvaluator/MathLinkTableEvaluatorBase.h>

void MathLinkTableEvaluatorBase::Replace_StarHat_With_e(std::string& str)
{
  static const std::string str_star_hat = "*^";
  static const std::string str_e        = "e";

  const std::string::size_type pos = str.find(str_star_hat);

  if(pos != std::string::npos)
  {
    str.replace(pos, str_star_hat.length(), str_e);
  }
}


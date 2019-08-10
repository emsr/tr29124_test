
#include <sstream>
#include <iomanip>

#include <interop/clr/e_float_clr.h>

e_float_clr::cli::e_float::e_float(System::String^ sys_str)
{
  std::string str(static_cast<std::string::size_type>(sys_str->Length), static_cast<char>('\0'));

  ::cli::array<wchar_t>^ wc_array = sys_str->ToCharArray();

  std::string::iterator it = str.begin();

  for each(wchar_t wc in wc_array)
  {
    *it = static_cast<char>(wc);
    ++it;
  }

  delete my_p;
  my_p = new ::e_float(str);
}

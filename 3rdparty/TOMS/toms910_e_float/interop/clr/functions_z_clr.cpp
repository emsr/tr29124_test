
#include <interop/clr/functions_z_clr.h>

e_float_clr::cli::ef_complex^ e_float_clr::cli::efz::chebyshev_t(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp)
{
  std::vector< ::ef_complex> vv;
  cli::ef_complex^ y = gcnew cli::ef_complex(::efz::chebyshev_t(n, z, &vv));

  vp->Clear();
  for(std::vector< ::ef_complex>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::ef_complex(*i));
  }

  return y;
}

e_float_clr::cli::ef_complex^ e_float_clr::cli::efz::chebyshev_u(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp)
{
  std::vector< ::ef_complex> vv;
  cli::ef_complex^ y = gcnew cli::ef_complex(::efz::chebyshev_u(n, z, &vv));

  vp->Clear();
  for(std::vector< ::ef_complex>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::ef_complex(*i));
  }

  return y;
}

e_float_clr::cli::ef_complex^ e_float_clr::cli::efz::hermite(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp)
{
  std::vector< ::ef_complex> vv;
  cli::ef_complex^ y = gcnew cli::ef_complex(::efz::hermite(n, z, &vv));

  vp->Clear();
  for(std::vector< ::ef_complex>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::ef_complex(*i));
  }

  return y;
}

e_float_clr::cli::ef_complex^ e_float_clr::cli::efz::laguerre(const UINT32 n, cli::ef_complex^ z, System::Collections::Generic::List<cli::ef_complex^>^ vp)
{
  std::vector< ::ef_complex> vv;
  cli::ef_complex^ y = gcnew cli::ef_complex(::efz::laguerre(n, z, &vv));

  vp->Clear();
  for(std::vector< ::ef_complex>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::ef_complex(*i));
  }

  return y;
}

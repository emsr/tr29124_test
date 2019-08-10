
#include <interop/clr/functions_clr.h>

System::Collections::Generic::List<e_float_clr::cli::e_float^>^ e_float_clr::cli::ef::airy_a_zero(const UINT32 k)
{
  System::Collections::Generic::List<cli::e_float^>^ zeros = gcnew System::Collections::Generic::List<cli::e_float^>;

  const std::deque< ::e_float> zd = ::ef::airy_a_zero(k);

  for(std::deque< ::e_float>::const_iterator i = zd.begin(); i != zd.end(); i++)
  {
    zeros->Add(gcnew cli::e_float(*i));
  }
  return zeros;
}

System::Collections::Generic::List<e_float_clr::cli::e_float^>^ e_float_clr::cli::ef::airy_b_zero(const UINT32 k)
{
  System::Collections::Generic::List<cli::e_float^>^ zeros = gcnew System::Collections::Generic::List<cli::e_float^>;

  const std::deque< ::e_float> zd = ::ef::airy_b_zero(k);

  for(std::deque< ::e_float>::const_iterator i = zd.begin(); i != zd.end(); i++)
  {
    zeros->Add(gcnew cli::e_float(*i));
  }

  return zeros;
}

System::Collections::Generic::List<e_float_clr::cli::e_float^>^ e_float_clr::cli::ef::cyl_bessel_j_zero(const INT32 n, const UINT32 k)
{
  System::Collections::Generic::List<cli::e_float^>^ zeros = gcnew System::Collections::Generic::List<cli::e_float^>;

  const std::deque< ::e_float> zd = ::ef::cyl_bessel_j_zero(n, k);

  for(std::deque< ::e_float>::const_iterator i = zd.begin(); i != zd.end(); i++)
  {
    zeros->Add(gcnew cli::e_float(*i));
  }

  return zeros;
}

System::Collections::Generic::List<e_float_clr::cli::e_float^>^ e_float_clr::cli::ef::cyl_bessel_j_zero(cli::e_float^ v, const UINT32 k)
{
  System::Collections::Generic::List<cli::e_float^>^ zeros = gcnew System::Collections::Generic::List<cli::e_float^>;

  const std::deque< ::e_float> zd = ::ef::cyl_bessel_j_zero(v, k);

  for(std::deque< ::e_float>::const_iterator i = zd.begin(); i != zd.end(); i++)
  {
    zeros->Add(gcnew cli::e_float(*i));
  }

  return zeros;
}

e_float_clr::cli::e_float^ e_float_clr::cli::ef::hyperg_pfq(System::Collections::Generic::List<cli::e_float^>^ a, System::Collections::Generic::List<cli::e_float^>^ b, cli::e_float^ x)
{
  std::deque< ::e_float> ad;
  std::deque< ::e_float> bd;

  for each(cli::e_float^ f in a)
  {
    ad.push_back(f->my_ref());
  }

  for each(cli::e_float^ f in b)
  {
    bd.push_back(f->my_ref());
  }

  return gcnew cli::e_float(::ef::hyperg_pfq(ad, bd, x));
}

System::Collections::Generic::List<System::UInt32>^ e_float_clr::cli::ef::prime(const UINT32 n)
{
  std::deque<UINT32> pd;
  ::ef::prime(n, pd);

  System::Collections::Generic::List<System::UInt32>^ p = gcnew System::Collections::Generic::List<System::UInt32>;
  for(std::deque<UINT32>::const_iterator i = pd.begin(); i != pd.end(); i++)
  {
    p->Add(System::UInt32(*i));
  }

  return p;
}

e_float_clr::cli::e_float^ e_float_clr::cli::ef::chebyshev_t(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp)
{
  std::vector< ::e_float> vv;
  cli::e_float^ y = gcnew cli::e_float(::ef::chebyshev_t(n, x, &vv));

  vp->Clear();
  for(std::vector< ::e_float>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::e_float(*i));
  }

  return y;
}

e_float_clr::cli::e_float^ e_float_clr::cli::ef::chebyshev_u(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp)
{
  std::vector< ::e_float> vv;
  cli::e_float^ y = gcnew cli::e_float(::ef::chebyshev_u(n, x, &vv));

  vp->Clear();
  for(std::vector< ::e_float>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::e_float(*i));
  }

  return y;
}

e_float_clr::cli::e_float^ e_float_clr::cli::ef::hermite(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp)
{
  std::vector< ::e_float> vv;
  cli::e_float^ y = gcnew cli::e_float(::ef::hermite(n, x, &vv));

  vp->Clear();
  for(std::vector< ::e_float>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::e_float(*i));
  }

  return y;
}

e_float_clr::cli::e_float^ e_float_clr::cli::ef::laguerre(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp)
{
  std::vector< ::e_float> vv;
  cli::e_float^ y = gcnew cli::e_float(::ef::laguerre(n, x, &vv));

  vp->Clear();
  for(std::vector< ::e_float>::const_iterator i = vv.begin(); i != vv.end(); i++)
  {
    vp->Add(gcnew cli::e_float(*i));
  }

  return y;
}

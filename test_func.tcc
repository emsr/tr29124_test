

template<typename Tp>
class filename_end
{
public:
  static const std::string suffix() { return std::string(".txt"); }
};

template<>
class filename_end<float>
{
public:
  static const std::string suffix() { return std::string("_f.txt"); }
};

template<>
class filename_end<double>
{
public:
  static const std::string suffix() { return std::string("_d.txt"); }
};

template<>
class filename_end<long double>
{
public:
  static const std::string suffix() { return std::string("_l.txt"); }
};


///
///  @brief  Fill an array with evenly spaces values between two limits.
///
template<typename Tp>
std::vector<Tp> fill_argument(const std::pair<Tp,Tp> & range,
                              const std::pair<bool,bool> & inclusive,
                              const unsigned int num_steps = 101)
{
  std::vector<Tp> argument;

  for (unsigned int i = 0; i < num_steps; ++i)
    {
      if (i == 0 && ! inclusive.first)
        {
          continue;
        }
      if (i == num_steps - 1 && ! inclusive.second)
        {
          continue;
        }

      Tp x = range.first + i * (range.second - range.first) / (num_steps - 1);

      argument.push_back(x);
    }

  return argument;
}


///
///  @brief  Run a function that takes one (1) argument.
///
template<typename Tp, typename Tp1>
void runtest(Tp function(const Tp1),
             const std::string & basename,
             const std::vector<Tp1> & argument1,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int i = 0; i < argument1.size(); ++i)
    {
      Tp1 x = argument1[i];

      output << "  " << std::setw(width) << x;

      try
        {
          Tp f = function(x);
          output << "  " << std::setw(width) << f;
        }
      catch (std::domain_error & err)
        {
          if (verbose)
            output << std::endl << err.what() << std::endl;
          continue;
        }
      catch (std::runtime_error & err)
        {
          if (verbose)
            output << std::endl << err.what() << std::endl;
          continue;
        }
      output << std::endl;
    }

  return;
}


///
///  @brief  Run a function that takes two (2) arguments.
///
template<typename Tp, typename Tp1, typename Tp2>
void runtest(Tp function(const Tp1, const Tp2),
             const std::string & basename,
             const std::vector<Tp1> & argument1,
             const std::vector<Tp2> & argument2,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int j = 0; j < argument2.size(); ++j)
    {
      Tp2 y = argument2[j];

      output << "  arg2 = " << std::setw(width) << y << std::endl;

      output << "  " << std::setw(width) << "arg1 = ";
      for (unsigned int i = 0; i < argument1.size(); ++i)
        {
          Tp1 x = argument1[i];

          output << "  " << std::setw(width) << x;
        }
      output << std::endl;

      for (unsigned int i = 0; i < argument1.size(); ++i)
        {
          Tp1 x = argument1[i];

          try
            {
              Tp f = function(x, y);
              output << "  " << std::setw(width) << f;
            }
          catch (std::domain_error & err)
            {
              if (verbose)
                output << std::endl << err.what() << std::endl;
              continue;
            }
          catch (std::runtime_error & err)
            {
              if (verbose)
                output << std::endl << err.what() << std::endl;
              continue;
            }
        }
      output << std::endl;
    }

  return;
}


///
///  @brief  Run a function that takes three (3) arguments.
///
template<typename Tp, typename Tp1, typename Tp2, typename Tp3>
void runtest(Tp function(const Tp1, const Tp2, const Tp3),
             const std::string & basename,
             const std::vector<Tp1> & argument1,
             const std::vector<Tp2> & argument2,
             const std::vector<Tp3> & argument3,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int k = 0; k < argument3.size(); ++k)
    {
      Tp3 z = argument3[k];
      output << "  arg3 = " << std::setw(width) << z << std::endl;

      for (unsigned int j = 0; j < argument2.size(); ++j)
        {
          Tp2 y = argument2[j];
          output << "  arg2 = " << std::setw(width) << y << std::endl;

          output << "  " << std::setw(width) << "arg1 = ";
          for (unsigned int i = 0; i < argument1.size(); ++i)
            {
              Tp1 x = argument1[i];

              output << "  " << std::setw(width) << x;
            }
          output << std::endl;

          for (unsigned int i = 0; i < argument1.size(); ++i)
            {
              Tp1 x = argument1[i];

              try
                {
                  Tp f = function(x, y, z);
                  output << "  " << std::setw(width) << f;
                }
              catch (std::domain_error & err)
                {
                  if (verbose)
                    output << std::endl << err.what() << std::endl;
                  continue;
                }
              catch (std::runtime_error & err)
                {
                  if (verbose)
                    output << std::endl << err.what() << std::endl;
                  continue;
                }
            }
          output << std::endl;
        }
      output << std::endl;
    }

  return;
}


///
///
///
template<typename Tp, typename Tp1, typename Tp2, typename Tp3, typename Tp4>
void runtest(Tp function(const Tp1, const Tp2, const Tp3, const Tp4),
             const std::string & basename,
             const std::vector<Tp1> & argument1,
             const std::vector<Tp2> & argument2,
             const std::vector<Tp3> & argument3,
             const std::vector<Tp4> & argument4,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int l = 0; l < argument4.size(); ++l)
    {
      Tp4 w = argument4[l];
      output << "  arg4 = " << std::setw(width) << w << std::endl;

      for (unsigned int k = 0; k < argument3.size(); ++k)
        {
          Tp3 z = argument3[k];
          output << "  arg3 = " << std::setw(width) << z << std::endl;

          for (unsigned int j = 0; j < argument2.size(); ++j)
            {
              Tp2 y = argument2[j];
              output << "  arg2 = " << std::setw(width) << y << std::endl;

              output << "  " << std::setw(width) << "arg1 = ";
              for (unsigned int i = 0; i < argument1.size(); ++i)
                {
                  Tp1 x = argument1[i];

                  output << "  " << std::setw(width) << x;
                }
              output << std::endl;

              for (unsigned int i = 0; i < argument1.size(); ++i)
                {
                  Tp1 x = argument1[i];

                  try
                    {
                      Tp f = function(x, y, z, w);
                      output << "  " << std::setw(width) << f;
                    }
                  catch (std::domain_error & err)
                    {
                      if (verbose)
                        output << std::endl << err.what() << std::endl;
                      continue;
                    }
                  catch (std::runtime_error & err)
                    {
                      if (verbose)
                        output << std::endl << err.what() << std::endl;
                      continue;
                    }
                }
              output << std::endl;
            }
          output << std::endl;
        }
      output << std::endl;
    }

  return;
}


///
///  @brief  Difference two one-argument functions.
///
template<typename Tp, typename Tp1>
void rundiff(Tp function1(const Tp1),
             Tp function2(const Tp1),
             const std::string & basename, const std::string & arg1,
             const std::vector<Tp1> & argument1,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  output << "  " << std::setw(width) << arg1;
  output << "  " << std::setw(width) << "f1";
  output << "  " << std::setw(width) << "f2";
  output << "  " << std::setw(width) << "delta";
  output << "  " << std::setw(width) << "frac";
  output << std::endl;

  Tp max_abs_diff = -Tp(1);
  Tp max_abs_frac = -Tp(1);
  for (unsigned int i = 0; i < argument1.size(); ++i)
    {
      Tp1 x = argument1[i];

      try
        {
          Tp f1 = function1(x);
          Tp f2 = function2(x);
          Tp diff = f1 - f2;
          if (std::abs(diff) > max_abs_diff)
            max_abs_diff = std::abs(diff);
          output << "  " << std::setw(width) << x;
          output << "  " << std::setw(width) << f1;
          output << "  " << std::setw(width) << f2;
          output << "  " << std::setw(width) << diff;
          if (std::abs(f2) > 0)
            {
              Tp frac = diff / f2;
              output << "  " << std::setw(width) << frac;
              if (std::abs(frac) > max_abs_frac)
                max_abs_frac = std::abs(frac);
            }
          else
            output << "  " << std::setw(width) << "-";
          output << std::endl;
        }
      catch (std::domain_error & err)
        {
          if (verbose)
            output << std::endl << err.what() << std::endl;
          continue;
        }
      catch (std::runtime_error & err)
        {
          if (verbose)
            output << std::endl << err.what() << std::endl;
          continue;
        }
    }
  if (max_abs_diff >= Tp(0))
    output << "max(abs(diff)) = " << max_abs_diff << std::endl;
  else
    output << "max(abs(diff)) = -" << std::endl;
  if (max_abs_frac >= Tp(0))
    output << "max(abs(frac)) = " << max_abs_frac << std::endl;
  else
    output << "max(abs(frac)) = -" << std::endl;

  return;
}


///
///  @brief  Difference two two-argument functions.
///
template<typename Tp, typename Tp1, typename Tp2>
void rundiff(Tp function1(const Tp1, const Tp2),
             Tp function2(const Tp1, const Tp2),
             const std::string & basename,
             const std::string & arg1, const std::vector<Tp1> & argument1,
             const std::string & arg2, const std::vector<Tp2> & argument2,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int i = 0; i < argument1.size(); ++i)
    {
      Tp1 x = argument1[i];
      output << "  " << arg1 << " = " << std::setw(width) << x << std::endl;

      output << "  " << std::setw(width) << arg2;
      output << "  " << std::setw(width) << "f1";
      output << "  " << std::setw(width) << "f2";
      output << "  " << std::setw(width) << "delta";
      output << "  " << std::setw(width) << "frac";
      output << std::endl;

      Tp max_abs_diff = -Tp(1);
      Tp max_abs_frac = -Tp(1);
      for (unsigned int j = 0; j < argument2.size(); ++j)
        {
          Tp2 y = argument2[j];

          try
            {
              Tp f1 = function1(x, y);
              Tp f2 = function2(x, y);
              Tp diff = f1 - f2;
              if (std::abs(diff) > max_abs_diff)
                max_abs_diff = std::abs(diff);
              output << "  " << std::setw(width) << y;
              output << "  " << std::setw(width) << f1;
              output << "  " << std::setw(width) << f2;
              output << "  " << std::setw(width) << diff;
              if (std::abs(f2) > 0)
                {
                  Tp frac = diff / f2;
                  output << "  " << std::setw(width) << frac;
                  if (std::abs(frac) > max_abs_frac)
                    max_abs_frac = std::abs(frac);
                }
              else
                output << "  " << std::setw(width) << "-";
              output << std::endl;
            }
          catch (std::domain_error & err)
            {
              if (verbose)
                output << std::endl << err.what() << std::endl;
              continue;
            }
          catch (std::runtime_error & err)
            {
              if (verbose)
                output << std::endl << err.what() << std::endl;
              continue;
            }
        }
      if (max_abs_diff >= Tp(0))
        output << "max(abs(diff)) = " << max_abs_diff << std::endl;
      else
        output << "max(abs(diff)) = -" << std::endl;
      if (max_abs_frac >= Tp(0))
        output << "max(abs(frac)) = " << max_abs_frac << std::endl;
      else
        output << "max(abs(frac)) = -" << std::endl;
      output << std::endl;
    }

  return;
}


///
///  @brief  Difference two three-argument functions.
///
template<typename Tp, typename Tp1, typename Tp2, typename Tp3>
void rundiff(Tp function1(const Tp1, const Tp2, const Tp3),
             Tp function2(const Tp1, const Tp2, const Tp3),
             const std::string & basename,
             const std::string & arg1, const std::vector<Tp1> & argument1,
             const std::string & arg2, const std::vector<Tp2> & argument2,
             const std::string & arg3, const std::vector<Tp3> & argument3,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int i = 0; i < argument1.size(); ++i)
    {
      Tp1 x = argument1[i];
      output << "  " << arg1 << " = " << std::setw(width) << x << std::endl;

      for (unsigned int j = 0; j < argument2.size(); ++j)
        {
          Tp2 y = argument2[j];
          output << "  " << arg2 << " = " << std::setw(width) << y << std::endl;

          output << "  " << std::setw(width) << arg3;
          output << "  " << std::setw(width) << "f1";
          output << "  " << std::setw(width) << "f2";
          output << "  " << std::setw(width) << "delta";
          output << "  " << std::setw(width) << "frac";
          output << std::endl;

          Tp max_abs_diff = -Tp(1);
          Tp max_abs_frac = -Tp(1);
          for (unsigned int k = 0; k < argument3.size(); ++k)
            {
              Tp3 z = argument3[k];

              try
                {
                  Tp f1 = function1(x, y, z);
                  Tp f2 = function2(x, y, z);
                  Tp diff = f1 - f2;
                  if (std::abs(diff) > max_abs_diff)
                    max_abs_diff = std::abs(diff);
                  output << "  " << std::setw(width) << z;
                  output << "  " << std::setw(width) << f1;
                  output << "  " << std::setw(width) << f2;
                  output << "  " << std::setw(width) << diff;
                  if (std::abs(f2) > 0)
                    {
                      Tp frac = diff / f2;
                      output << "  " << std::setw(width) << frac;
                      if (std::abs(frac) > max_abs_frac)
                        max_abs_frac = std::abs(frac);
                    }
                  else
                    output << "  " << std::setw(width) << "-";
                  output << std::endl;
                }
              catch (std::domain_error & err)
                {
                  if (verbose)
                    output << std::endl << err.what() << std::endl;
                  continue;
                }
              catch (std::runtime_error & err)
                {
                  if (verbose)
                    output << std::endl << err.what() << std::endl;
                  continue;
                }
            }
          if (max_abs_diff >= Tp(0))
            output << "max(abs(diff)) = " << max_abs_diff << std::endl;
          else
            output << "max(abs(diff)) = -" << std::endl;
          if (max_abs_frac >= Tp(0))
            output << "max(abs(frac)) = " << max_abs_frac << std::endl;
          else
            output << "max(abs(frac)) = -" << std::endl;
          output << std::endl;
        }
      output << std::endl;
    }

  return;
}


///
///  @brief  Difference two four-argument functions.
///
template<typename Tp, typename Tp1, typename Tp2, typename Tp3, typename Tp4>
void rundiff(Tp function1(const Tp1, const Tp2, const Tp3, const Tp4),
             Tp function2(const Tp1, const Tp2, const Tp3, const Tp4),
             const std::string & basename,
             const std::string & arg1, const std::vector<Tp1> & argument1,
             const std::string & arg2, const std::vector<Tp2> & argument2,
             const std::string & arg3, const std::vector<Tp3> & argument3,
             const std::string & arg4, const std::vector<Tp4> & argument4,
             const bool verbose = false)
{
  std::string filename = basename + filename_end<Tp>::suffix();

  std::ofstream output(filename.c_str());
  output.precision(std::numeric_limits<Tp>::digits10);
  output.flags(std::ios::showpoint);
  int width = output.precision() + 6;

  for (unsigned int i = 0; i < argument1.size(); ++i)
    {
      Tp1 w = argument1[i];
      output << "  " << arg1 << " = " << std::setw(width) << w << std::endl;

      for (unsigned int j = 0; j < argument2.size(); ++j)
        {
          Tp2 x = argument2[j];
          output << "  " << arg2 << " = " << std::setw(width) << x << std::endl;

          for (unsigned int k = 0; k < argument3.size(); ++k)
            {
              Tp3 y = argument3[k];
              output << "  " << arg3 << " = " << std::setw(width) << y << std::endl;

              output << "  " << std::setw(width) << arg4;
              output << "  " << std::setw(width) << "f1";
              output << "  " << std::setw(width) << "f2";
              output << "  " << std::setw(width) << "delta";
              output << "  " << std::setw(width) << "frac";
              output << std::endl;

              Tp max_abs_diff = -Tp(1);
              Tp max_abs_frac = -Tp(1);
              for (unsigned int l = 0; l < argument4.size(); ++l)
                {
                  Tp4 z = argument4[l];

                  try
                    {
                      Tp f1 = function1(w, x, y, z);
                      Tp f2 = function2(w, x, y, z);
                      Tp diff = f1 - f2;
                      if (std::abs(diff) > max_abs_diff)
                        max_abs_diff = std::abs(diff);
                      output << "  " << std::setw(width) << z;
                      output << "  " << std::setw(width) << f1;
                      output << "  " << std::setw(width) << f2;
                      output << "  " << std::setw(width) << diff;
                      if (std::abs(f2) > 0)
                        {
                          Tp frac = diff / f2;
                          output << "  " << std::setw(width) << frac;
                          if (std::abs(frac) > max_abs_frac)
                            max_abs_frac = std::abs(frac);
                        }
                      else
                        output << "  " << std::setw(width) << "-";
                      output << std::endl;
                    }
                  catch (std::domain_error & err)
                    {
                      if (verbose)
                        output << std::endl << err.what() << std::endl;
                      continue;
                    }
                  catch (std::runtime_error & err)
                    {
                      if (verbose)
                        output << std::endl << err.what() << std::endl;
                      continue;
                    }
                }
              if (max_abs_diff >= Tp(0))
                output << "max(abs(diff)) = " << max_abs_diff << std::endl;
              else
                output << "max(abs(diff)) = -" << std::endl;
              if (max_abs_frac >= Tp(0))
                output << "max(abs(frac)) = " << max_abs_frac << std::endl;
              else
                output << "max(abs(frac)) = -" << std::endl;
              output << std::endl;
            }
          output << std::endl;
        }
      output << std::endl;
    }

  return;
}


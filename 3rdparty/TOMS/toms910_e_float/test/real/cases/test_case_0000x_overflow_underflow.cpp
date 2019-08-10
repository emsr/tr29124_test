
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCaseOverflowUnderflowBase : public TestCaseReal
    {
    protected:

      mutable bool my_test_result;

      TestCaseOverflowUnderflowBase() : my_test_result(false) { }

    public:

      virtual ~TestCaseOverflowUnderflowBase() { }

      virtual bool execute(const bool b_write_output) const
      {
        std::cout << name() << " : ";

        std::vector<e_float> e_float_data;

        // Calculate the e_float test data.
        e_float_test(e_float_data);

        // Optionally write the e_float test data to an output file.
        if(b_write_output)
        {
          if(!write_output_file(e_float_data))
          {
            std::cout << "Can not write output: FAIL" << std::endl;
            return false;
          }
        }

        if(my_test_result)
        {
          std::cout << "Numerical compare OK: PASS"  << std::endl;
          return true;
        }
        else
        {
          std::cout << "Numerical compare not OK: FAIL"  << std::endl;
          return false;
        }
      }
    };

    class TestCase_case_00001_overflow_mul_x : public TestCaseOverflowUnderflowBase
    {
    public:
      TestCase_case_00001_overflow_mul_x() { }
      virtual ~TestCase_case_00001_overflow_mul_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00001_overflow_mul_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();

        e_float y = ef::googol() * ef::euler_gamma();

        data.push_back(y);

        INT32 k;
        for(k = static_cast<INT32>(0); k < static_cast<INT32>(1000); k++)
        {
          y = y * y;

          data.push_back(y);

          if(ef::isinf(y))
          {
            break;
          }
        }

        my_test_result = (k > static_cast<INT32>(1)) && (k < static_cast<INT32>(1000));
      }
    };

    class TestCase_case_00002_underflow_mul_x : public TestCaseOverflowUnderflowBase
    {
    public:
      TestCase_case_00002_underflow_mul_x() { }
      virtual ~TestCase_case_00002_underflow_mul_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00002_underflow_mul_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();

        e_float y = ef::one() / (ef::googol() * ef::euler_gamma());

        data.push_back(y);

        INT32 k;
        for(k = static_cast<INT32>(0); k < static_cast<INT32>(1000); k++)
        {
          y = y * y;

          data.push_back(y);

          if(ef::iszero(y))
          {
            break;
          }
        }

        my_test_result = (k > static_cast<INT32>(1)) && (k < static_cast<INT32>(1000));
      }
    };

    class TestCase_case_00003_overflow_x_mul_by_n : public TestCaseOverflowUnderflowBase
    {
    public:
      TestCase_case_00003_overflow_x_mul_by_n() { }
      virtual ~TestCase_case_00003_overflow_x_mul_by_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00003_overflow_x_mul_by_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();

        e_float y = std::numeric_limits<e_float>::max() / static_cast<INT32>(100);

        data.push_back(y);

        INT32 k;
        for(k = static_cast<INT32>(0); k < static_cast<INT32>(1000); k++)
        {
          y = y * static_cast<INT32>(2);

          data.push_back(y);

          if(ef::isinf(y))
          {
            break;
          }
        }

        my_test_result = (k > static_cast<INT32>(1)) && (k < static_cast<INT32>(1000));
      }
    };

    class TestCase_case_00004_underflow_x_div_by_n : public TestCaseOverflowUnderflowBase
    {
    public:
      TestCase_case_00004_underflow_x_div_by_n() { }
      virtual ~TestCase_case_00004_underflow_x_div_by_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00004_underflow_x_div_by_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();

        e_float y = std::numeric_limits<e_float>::min() * static_cast<INT32>(100);

        data.push_back(y);

        INT32 k;
        for(k = static_cast<INT32>(0); k < static_cast<INT32>(1000); k++)
        {
          y = y / static_cast<INT32>(2);

          data.push_back(y);

          if(ef::iszero(y))
          {
            break;
          }
        }

        my_test_result = (k > static_cast<INT32>(1)) && (k < static_cast<INT32>(1000));
      }
    };

    bool test_case_00001_overflow_mul_x (const bool b_write_output)
    {
      return TestCase_case_00001_overflow_mul_x().execute(b_write_output);
    }
    bool test_case_00002_underflow_mul_x(const bool b_write_output)
    {
      return TestCase_case_00002_underflow_mul_x().execute(b_write_output);
    }
    bool test_case_00003_overflow_x_mul_by_n(const bool b_write_output)
    {
      return TestCase_case_00003_overflow_x_mul_by_n().execute(b_write_output);
    }
    bool test_case_00004_underflow_x_div_by_n(const bool b_write_output)
    {
      return TestCase_case_00004_underflow_x_div_by_n().execute(b_write_output);
    }
  }
}

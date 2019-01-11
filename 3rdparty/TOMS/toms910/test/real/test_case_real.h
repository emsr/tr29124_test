#ifndef _TEST_CASE_REAL_2009_10_24_H_
  #define _TEST_CASE_REAL_2009_10_24_H_

  #include <test/test_case_base.h>

  namespace test
  {
    namespace real
    {
      class TestCaseReal : public TestCaseBase<e_float>
      {
      private:

        TestCaseReal(const TestCaseReal&);
        const TestCaseReal& operator=(const TestCaseReal&);

      protected:

        TestCaseReal() { }

      public:

        virtual ~TestCaseReal() { }
      };
    }
  }

#endif // _TEST_CASE_REAL_2009_10_24_H_

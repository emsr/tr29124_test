#ifndef _TEST_CASE_IMAG_2009_10_24_H_
  #define _TEST_CASE_IMAG_2009_10_24_H_

  #include <test/test_case_base.h>

  namespace test
  {
    namespace imag
    {
      class TestCaseImag : public TestCaseBase<ef_complex>
      {
      private:

        TestCaseImag(const TestCaseImag&);
        const TestCaseImag& operator=(const TestCaseImag&);

      protected:

        TestCaseImag() { }

      public:

        virtual ~TestCaseImag() { }
      };
    }
  }

#endif // _TEST_CASE_IMAG_2009_10_24_H_

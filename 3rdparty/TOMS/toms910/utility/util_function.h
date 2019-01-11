
#ifndef _UTIL_FUNCTION_2009_10_27_H_
  #define _UTIL_FUNCTION_2009_10_27_H_

  namespace Util
  {
    template<typename T> class Function
    {
    private:

      const Function& operator=(const Function&);
      Function(const Function&);

    protected:

      Function() { }

    public:

      virtual ~Function() { }

    private:

      virtual T my_function(const T& x) const = 0;

    protected:

      T function(const T& x) const { return my_function(x); }
    };
  }

#endif // _UTIL_FUNCTION_2009_10_27_H_

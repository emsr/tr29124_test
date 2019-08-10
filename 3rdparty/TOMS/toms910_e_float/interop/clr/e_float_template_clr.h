
#ifndef _E_FLOAT_TEMPLATE_CLR_2010_01_06_H_
  #define _E_FLOAT_TEMPLATE_CLR_2010_01_06_H_

  namespace e_float_clr
  {
    namespace cli
    {
      template<typename T> public ref class e_float_template
      {
      protected:

        T* my_p;

        e_float_template() : my_p(new T()) { }

        explicit e_float_template(const  INT32 n) : my_p(new T(n)) { }
        explicit e_float_template(const  INT64 n) : my_p(new T(n)) { }
        explicit e_float_template(const UINT32 u) : my_p(new T(u)) { }
        explicit e_float_template(const UINT64 u) : my_p(new T(u)) { }
        explicit e_float_template(const double d) : my_p(new T(d)) { }

        explicit e_float_template(const T& t) : my_p(new T(t)) { }
        e_float_template(e_float_template^ t) : my_p(new T(*(t->my_p))) { }

      public:

        virtual ~e_float_template() { this->!e_float_template(); }

        const T& my_ref(void) { return *my_p; }
        const T* my_ptr(void) { return  my_p; }

        operator const T&() { return my_ref(); }

      private:

        !e_float_template() { delete my_p; my_p = static_cast<T*>(0u); }

      public:

        e_float_template^ operator =(e_float_template^ v) { (*my_p)  = (*(v->my_p)); return this; }
        e_float_template^ operator+=(e_float_template^ v) { (*my_p) += (*(v->my_p)); return this; }
        e_float_template^ operator-=(e_float_template^ v) { (*my_p) -= (*(v->my_p)); return this; }
        e_float_template^ operator*=(e_float_template^ v) { (*my_p) *= (*(v->my_p)); return this; }
        e_float_template^ operator/=(e_float_template^ v) { (*my_p) /= (*(v->my_p)); return this; }

        e_float_template^ operator+=(const INT32 n) { (*my_p) += n; return this; }
        e_float_template^ operator-=(const INT32 n) { (*my_p) -= n; return this; }
        e_float_template^ operator*=(const INT32 n) { (*my_p) *= n; return this; }
        e_float_template^ operator/=(const INT32 n) { (*my_p) /= n; return this; }

        e_float_template^ negate(void) { (*my_p) = -(*my_p); return this; }

        e_float_template^ operator++(void) { ++(*my_p); return this; }
        e_float_template^ operator--(void) { --(*my_p); return this; }
      };
    }
  }

#endif // _E_FLOAT_TEMPLATE_CLR_2010_01_06_H_

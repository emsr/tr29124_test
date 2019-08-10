
#ifndef _E_FLOAT_PYD_ADAPT_2010_01_05_H_
  #define _E_FLOAT_PYD_ADAPT_2010_01_05_H_

  #include <boost/python/stl_iterator.hpp>

  namespace e_float_pyd
  {
    namespace pyd
    {
      namespace adapt
      {
        template<typename T> inline boost::python::list std_deque_to_pylist(const std::deque<T>& std_deq)
        {
          boost::python::list pylist;
          for(typename std::deque<T>::const_iterator it = std_deq.begin(); it != std_deq.end(); it++)
          {
            pylist.append(*it);
          }
          return pylist;
        }

        template<typename T> inline boost::python::list std_vector_to_pylist(const std::vector<T>& std_vec)
        {
          boost::python::list pylist;
          for(typename std::vector<T>::const_iterator it = std_vec.begin(); it != std_vec.end(); it++)
          {
            pylist.append(*it);
          }
          return pylist;
        }

        template<typename T> inline std::deque<T> pylist_to_std_deque(const boost::python::list& pylist)
        {
          const typename boost::python::stl_input_iterator<T> begin(pylist), end;
          return std::deque<T>(begin, end);
        }

        template<typename T> inline std::vector<T> pylist_to_std_vector(const boost::python::list& pylist)
        {
          const typename boost::python::stl_input_iterator<T> begin(pylist), end;
          return std::vector<T>(begin, end);
        }
      }
    }
  }

#endif // _E_FLOAT_PYD_ADAPT_2010_01_05_H_

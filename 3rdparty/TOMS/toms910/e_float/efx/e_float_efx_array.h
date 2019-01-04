
#ifndef _E_FLOAT_EFX_ARRAY_2008_05_25_H_
  #define _E_FLOAT_EFX_ARRAY_2008_05_25_H_

  #include <algorithm>
  #include <stdexcept>

  // Implement a highly speed-optimized array<T, N> class for the internal data-type
  // of efx::e_float. This implementation has been taken and adapted from boost::array<>.

  namespace efx
  {
    template<typename T, const std::size_t N> class array
    {
    public:

      T elems[N];    // fixed-size array of elements of type T

      // type definitions
      typedef T              value_type;
      typedef T*             iterator;
      typedef const T*       const_iterator;
      typedef T&             reference;
      typedef const T&       const_reference;
      typedef std::size_t    size_type;
      typedef std::ptrdiff_t difference_type;

      // reverse iterator support
      typedef std::reverse_iterator<iterator> reverse_iterator;
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

      // iterator support
      iterator begin(void)             { return elems; }
      const_iterator begin(void) const { return elems; }
      iterator end(void)               { return elems + N; }
      const_iterator end(void) const   { return elems + N; }

      reverse_iterator rbegin(void)             { return reverse_iterator(elems + N); }
      const_reverse_iterator rbegin(void) const { return const_reverse_iterator(elems + N); }

      reverse_iterator rend(void)             { return reverse_iterator(elems); }
      const_reverse_iterator rend(void) const { return const_reverse_iterator(elems); }

      // operator[] without range check
      reference operator[](const size_type i)             { return elems[i]; }
      const_reference operator[](const size_type i) const { return elems[i]; }

      // at() with range check
      reference at(const size_type i)             { rangecheck(i); return elems[i]; }
      const_reference at(const size_type i) const { rangecheck(i); return elems[i]; }

      // front() and back()
      reference front(void)             { return elems[0u]; }
      const_reference front(void) const { return elems[0u]; }
      
      reference back(void)             { return elems[N - 1u]; }
      const_reference back(void) const { return elems[N - 1u]; }

      // size is constant
      static size_type size(void)     { return N; }
      static bool empty(void)         { return false; }
      static size_type max_size(void) { return N; }

      static const size_type static_size = N;

      // swap (note: linear complexity)
      void swap(array<T, N>& y) { static_cast<void>(std::swap_ranges(begin(), end(), y.begin())); }

      // direct access to data (read-only)
      const T* data(void) const { return elems; }
      T* data(void)             { return elems; }

      // use array as C array (direct read/write access to data)
      T* c_array(void) { return elems; }

      // assignment with type conversion
      template<typename T2> array<T, N>& operator=(const array<T2, N>& rhs)
      {
        static_cast<void>(std::copy(rhs.begin(), rhs.end(), begin()));
        return *this;
      }

      // assign one value to all elements
      void assign(const T& value)
      {
        std::fill_n(elems, N, value);
      }

      // Kormanyos: added comparison operators for array.
      bool operator==(const array& rhs) const { return std::equal(begin(), end(), rhs.begin()); }
      bool operator!=(const array& rhs) const { return !(operator==(rhs)); }
      bool operator< (const array& rhs) const { return std::lexicographical_compare(begin(), end(), rhs.begin(), rhs.end()); }
      bool operator> (const array& rhs) const { return  (rhs < *this); }
      bool operator<=(const array& rhs) const { return !(rhs < *this); }
      bool operator>=(const array& rhs) const { return !(*this < rhs); }

    private:

      // check range (may be private because it is static)
      static void rangecheck(const size_type i)
      {
        if(i >= N)
        {
          throw std::out_of_range("array<>: index out of range");
        }
      }

    };
  }

#endif // _E_FLOAT_EFX_ARRAY_2008_05_25_H_

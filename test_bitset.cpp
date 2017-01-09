
#include <bitset>

/*
Poor bitset needs:


      /**
       *  Use a subset of a string.
       *  @param  __s  A string of @a 0 and @a 1 characters.
       *  @param  __position  Index of the first character in @a __s to use;
       *                    defaults to zero.
       *  @throw  std::out_of_range  If @a pos is bigger the size of @a __s.
       *  @throw  std::invalid_argument  If a character appears in the string
       *                                 which is neither @a 0 nor @a 1.
       * /
      template<class _CharT, class _Traits>
	explicit
	bitset(std::basic_string_view<_CharT, _Traits> __s,
	       size_t __position = 0)
	: _Base()
	{
	  _M_check_initial_position(__s, __position);
	  _M_copy_from_string(__s, __position,
			      std::basic_string_view<_CharT, _Traits>::npos,
			      _CharT('0'), _CharT('1'));
	}

      // _GLIBCXX_RESOLVE_LIB_DEFECTS
      // 396. what are characters zero and one.
      template<class _CharT, class _Traits>
	bitset(std::basic_string_view<_CharT, _Traits> __s,
	       size_t __position, size_t __n,
	       _CharT __zero, _CharT __one = _CharT('1'))
	: _Base()
	{
	  _M_check_initial_position(__s, __position);
	  _M_copy_from_string(__s, __position, __n, __zero, __one);
	}

  template<char... _Bits>
    auto
    operator""bitset()
    -> std::bitset<sizeof...(_Bits)>
    {
      using __bset_t = std::bitset<sizeof...(_Bits)>;
      const char __bits[]{_Bits...,'\0'};
      return __bset_t(__bits);
    }
*/

template<char... _Bits>
  auto
  operator""_bits()
  -> std::bitset<sizeof...(_Bits)>
  {
    using __bset_t = std::bitset<sizeof...(_Bits)>;
    const char __bits[]{_Bits...,'\0'};
    return __bset_t(__bits);
  }

/*
  auto
  operator""_bits(const char* __bits, std::size_t __len)
  -> std::bitset<__len>
  {
    using __bset_t = std::bitset<__len>;
    return __bset_t(__bits);
  }
*/

int
main()
{
//11.
  auto bigpi = 0010010000111111011010101000100010000101101000110000100011010011_bits;
/*    "0010010000111111011010101000100010000101101000110000100011010011"
    "0001001100011001100010100010111000000011011100000111001101000100"
    "1010010000001001001110000010001000101001100111110011000111010000"
    "0000100000101110111110101001100011101100010011100110110010001001"
    "0100010100101000001000011110011000111000110100000001001101110111"
    "1011111001010100011001101100111100110100111010010000110001101100"
    "1100000010101100001010011011011111001001011111000101000011011101"
    "0011111110000100110101011011010110110101010001110000100100010111"
    "1001001000010110110101011101100110001001011110011111101100011011"
    "1101000100110001000010111010011010011000110111111011010110101100"
    "0010111111111101011100101101101111010000000110101101111110110111"
    "1011100011100001101011111110110101101010001001100111111010010110"
    "1011101001111100100100000100010111110001001011000111111110011001"
    "0010010010100001100110010100011110110011100100010110110011110111"
    "0000100000000001111100101110001010000101100011101111110000010110"
    "0110001101101001001000001101100001110001010101110100111001101001"_bits;*/
}

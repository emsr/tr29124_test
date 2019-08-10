// Synopsis of the e_float global arithmetic interface.
// Global operators post-increment and post-decrement
e_float operator++(e_float& u, int);
e_float operator--(e_float& u, int);

// Global unary operators of e_float reference.
      e_float  operator-(const e_float& u);
      e_float& operator+(      e_float& u);
const e_float& operator+(const e_float& u);

// Global add/sub/mul/div of const e_float reference with const e_float reference
e_float operator+(const e_float& u, const e_float& v);
e_float operator-(const e_float& u, const e_float& v);
e_float operator*(const e_float& u, const e_float& v);
e_float operator/(const e_float& u, const e_float& v);

// Specialization for global add/sub/mul/div of const e_float reference with INT32
e_float operator+(const e_float& u, const INT32 n);
e_float operator-(const e_float& u, const INT32 n);
e_float operator*(const e_float& u, const INT32 n);
e_float operator/(const e_float& u, const INT32 n);

e_float operator+(const INT32 n, const e_float& u);
e_float operator-(const INT32 n, const e_float& u);
e_float operator*(const INT32 n, const e_float& u);
e_float operator/(const INT32 n, const e_float& u);

// Specializations of global self-add/sub/mul-div of e_float reference with INT32
e_float& operator+=(e_float& u, const INT32 n);
e_float& operator-=(e_float& u, const INT32 n);
e_float& operator*=(e_float& u, const INT32 n);
e_float& operator/=(e_float& u, const INT32 n);

// Global comparison operators of const e_float reference with const e_float reference
bool operator< (const e_float& u, const e_float& v);
bool operator<=(const e_float& u, const e_float& v);
bool operator==(const e_float& u, const e_float& v);
bool operator!=(const e_float& u, const e_float& v);
bool operator>=(const e_float& u, const e_float& v);
bool operator> (const e_float& u, const e_float& v);

// Specializations of global comparison of e_float reference with INT32
bool operator< (const e_float& u, const INT32 n);
bool operator<=(const e_float& u, const INT32 n);
bool operator==(const e_float& u, const INT32 n);
bool operator!=(const e_float& u, const INT32 n);
bool operator>=(const e_float& u, const INT32 n);
bool operator> (const e_float& u, const INT32 n);

bool operator< (const INT32 n, const e_float& u);
bool operator<=(const INT32 n, const e_float& u);
bool operator==(const INT32 n, const e_float& u);
bool operator!=(const INT32 n, const e_float& u);
bool operator>=(const INT32 n, const e_float& u);
bool operator> (const INT32 n, const e_float& u);

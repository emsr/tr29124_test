
  if (aneg)
    {
      if (sint)
	{ /* ... s is an integer ... */
	  if (ztiny)
	    { /* ... Return first term of series. */
	      if (int(s) % 2 == 0)
		sign = 1;
	      else
		sign = -1;
	      result = sign * 1.0 / std::pow(std::abs(a), s);
	    }
	  else
	    { /* ... Transform a to positive. */
	      m = -int(std::floor(a));
	      a1 += m;
	      sum1 = 0.0;
	      if ((int) s % 2 == 0)
		sign = 1;
	      else
		sign = -1;
	      for (i = 0; i <= m-1; i++)
		{
		  sum1 += sign * std::pow(std::abs(z), i)
			       / std::pow(std::abs(a + i), s);
		  if (z < 0.0)
		    sign = -sign;
		}
	    }
	}
      else
	{ /* ... s is not an integer. (Return error because pow() is not defined.) */
	  result = 1.0;
	  iter = 0;
	  return 3;
	}
    }
  else
    {
      if (ztiny)
	{ /* ... a > 0. Return first term of series.) */
	  result = 1.0 / std::pow(a, s);
	  iter = 1;
	  return 0;
	}
    }

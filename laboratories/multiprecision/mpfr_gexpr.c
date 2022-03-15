/*****************************************************************************/
/* Extended gexpr.c for MPFR                                                 */
/* Original gexpr.c from http://swox.com/gmp/gexpr.c                         */
/* Version 0.1 by Tomonori Kouya, Jan. 25, 2004: "Nearly Correctly Rounded"  */
/* support suggested by Paul Zimmermann (Thanks!)                            */
/* Version 0.2 modified by Paul Zimmermann, Feb. 22, 2005.                   */
/* Version 0.3 modified by Tomonori Kouya, Sep. 22, 2005.                    */
/* Version 0.3a modified by Tomonori Kouya, Sep. 22, 2011                    */
/* Version 0.3b modified by Tomonori Kouya, Apr. 22, 2014                    */
/* ************************** CAUTION! ***************************************/
/* GNU MP and MPFR are based on LGPL (Lesser GNU Public License),            */
/* also BNCpack and MPIBNCpack. BUT the "mpfr_gexpr.c"'s copyright is GPL    */
/* because of GPLed original "gexpr.c."                                      */
/*                                                                           */
/* Expression evaluation using plain floating-point arithmetic.              */
/* Copyright (C) 1999, 2000, 2001, 2002, 2003 Free Software Foundation, Inc. */
/* This program is free software; you can redistribute it and/or modify it   */
/* under the terms of the GNU General Public License as published by the     */
/* Free Software Foundation; either version 2, or (at your option) any later */
/* version. This program is distributed in the hope that it will be useful,  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General */
/* Public License for more details.                                          */
/*                                                                           */
/* You should have received a copy of the GNU General Public License along   */
/* with this program; see the file COPYING.  If not, write to the Free       */
/* Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,  */
/* USA.                                                                      */
/*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h> /* for atoi, exit */
#include <setjmp.h>
#include <math.h>
#include <ctype.h> /* for isalnum */

#include <gmp.h>
#include <mpfr.h>

#if (MPFR_VERSION <= MPFR_VERSION_NUM(2,1,0))
#error "mpfr version should be >= 2.1.0"
#endif

jmp_buf top;

static char *bp, *expr_start_p, *expr_end_p;
static int ch;
static int previous_result_valid_flag;

mp_rnd_t default_rmode = GMP_RNDN; /* round to nearest */
mpfr_t previous_result;

static void mpfr_factor (mpfr_t);
static void mpfr_number (mpfr_t);

static void
skip (void)
{
  do
    ch = *bp++;
  while (ch == ' ' || ch == '\n');
}

static void
new_error (char *funcname)
{
  fprintf (stderr, "Syntax error in %s\n", funcname);
  longjmp (top, 1);
}

/* ret^expo() = exp(expo() * log(ret)) */
static void
mpfr_expo (mpfr_t ret)
{
  mpfr_t tmp;

  mpfr_init2 (tmp,  mpfr_get_prec (ret));

  mpfr_factor (ret);
  if (ch == '^')
    {
      skip ();
      mpfr_expo (tmp);

      // suggested by P.Zemmermann (Thanks!)
      mpfr_pow(ret, ret, tmp, default_rmode);
    //  mpfr_log (ret, ret, default_rmode); /* log (ret) */
    //  mpfr_mul (ret, tmp, ret, default_rmode);
    //  mpfr_exp (ret, ret, default_rmode); /* e^(expo() * log_e(ret))) */
    }

   mpfr_clear (tmp);
}

static void
mpfr_term (mpfr_t ret)
{
  mpfr_t tmp;

  mpfr_init2 (tmp, mpfr_get_prec(ret));

  mpfr_expo (ret);
  for (;;)
    switch (ch)
      {
      case '*':
	skip ();
	mpfr_expo (tmp);
	mpfr_mul (ret, ret, tmp, default_rmode);
	break;
      case '/':
	skip ();
	mpfr_expo (tmp);
	mpfr_div (ret, ret, tmp, default_rmode);
	break;
/*      case '%':
	skip ();
	t = fmod (t, expo ());
	break;
*/
      default:
  	 mpfr_clear (tmp);
         return;
      }
}

static void
mpfr_expr (mpfr_t ret)
{
  mpfr_t tmp;

  mpfr_init2 (tmp, mpfr_get_prec (ret));

  if (ch == '+')
    {
      skip ();
      mpfr_term (ret);
    }
  else if (ch == '-')
    {
      skip ();
      mpfr_term (ret);
      mpfr_neg (ret, ret, default_rmode);
    }
  else
    mpfr_term (ret);

  while (ch == '+' || ch == '-')
    {
      char op;
      op = ch;
      skip ();
      if (op == '-')
        {
          mpfr_term (tmp);
          mpfr_sub (ret, ret, tmp, default_rmode);
        }
      else
        {
          mpfr_term (tmp);
          mpfr_add (ret, ret, tmp, default_rmode);
        }
    }

  mpfr_clear (tmp);
}

struct log_functions
{
  char *spelling;
  int (* evalfn) (mpfr_ptr, int *, mpfr_srcptr, mpfr_rnd_t);
};

struct log_functions log_fns[] =
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,3,0))
  {"lgamma", mpfr_lgamma}
#endif
};

struct int_functions
{
  char *spelling;
  int (* evalfn) (mpfr_ptr, mpfr_srcptr);
};

struct int_functions int_fns[] =
{
  {"ceil", mpfr_ceil},
  {"floor", mpfr_floor}
};

struct int_arg_function
{
  char *spelling;
  int (* evalfn) (mpfr_ptr, long int, mpfr_srcptr, mp_rnd_t);
};

struct int_arg_function int_arg_fns[] =
{
  {"jn", mpfr_jn},
  {"yn", mpfr_yn}
};

struct twoarg_functions
{
  char *spelling;
  int (* evalfn) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t);
};

struct twoarg_functions twoarg_fns[] =
{
  {"agm", mpfr_agm},
  {"hypot", mpfr_hypot},
  {"atan2", mpfr_atan2},
  {"pow", mpfr_pow}
};

struct functions
{
  char *spelling;
  int (* evalfn) (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
};

struct functions fns[] =
{
  {"abs", mpfr_abs},
  {"acos", mpfr_acos},
  {"acosh", mpfr_acosh},
  {"asin", mpfr_asin},
  {"asinh", mpfr_asinh},
  {"atan", mpfr_atan},
  {"atanh", mpfr_atanh},
  {"cos", mpfr_cos},
  {"cosh", mpfr_cosh},
  {"exp", mpfr_exp},
  {"erf", mpfr_erf},
  {"expm1", mpfr_expm1},
  {"gamma", mpfr_gamma},
  {"lngamma", mpfr_lngamma},
  {"log", mpfr_log},
  {"log10", mpfr_log10},
  {"log2", mpfr_log2},
  {"sin", mpfr_sin},
  {"sinh", mpfr_sinh},
  {"sqrt", mpfr_sqrt},
  {"tan", mpfr_tan},
  {"tanh", mpfr_tanh},
  {"zeta", mpfr_zeta},
// 2005.0922
// appended functions to mpfr-2.2.0
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,2,0))
  {"cot", mpfr_cot},
  {"coth", mpfr_coth},
  {"csc", mpfr_csc},
  {"csch", mpfr_csch},
  {"eint", mpfr_eint},
  {"erfc", mpfr_erfc},
  {"sec", mpfr_sec},
  {"sech", mpfr_sech},
#endif
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,3,0))
  {"j0", mpfr_j0},
  {"j1", mpfr_j1},
  {"y0", mpfr_y0},
  {"y1", mpfr_y1},
#endif
  {"expm1", mpfr_expm1},
  {"exp10", mpfr_exp10},
  {"exp2", mpfr_exp2},
  {"li2", mpfr_li2}, // 2.4
#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
  {"ai", mpfr_ai}, // 3.0
  //{"bi", mpfr_bi},
  {"digamma", mpfr_digamma}, // 3.0
#endif
  {0, 0}
};


static void
mpfr_factor (mpfr_t ret)
{
  mpfr_t tmp, tmp2;

  mpfr_init2 (tmp, mpfr_get_prec (ret));
  mpfr_init2 (tmp2, mpfr_get_prec (ret));

  for (int i = 0; fns[i].spelling != 0; i++)
    {
      char *spelling = fns[i].spelling;
      int len = strlen (spelling);
      if (strncmp (spelling, bp - 1, len) == 0 && ! isalnum (bp[-1 + len]))
	{
	  bp += len - 1;

	  skip ();
	  if (ch != '(')
	    new_error ("mpfr_factor 1");

	  skip ();
	  mpfr_expr (tmp);

	  if (ch != ')')
	    new_error ("mpfr_factor 2");

	  skip ();

	  (fns[i].evalfn) (ret, tmp, default_rmode);

          mpfr_clear (tmp);

	  return;
	}
    }

  for (int i = 0; int_fns[i].spelling != 0; i++)
    {
      char *spelling = int_fns[i].spelling;
      int len = strlen (spelling);
      if (strncmp (spelling, bp - 1, len) == 0 && ! isalnum (bp[-1 + len]))
	{
	  bp += len - 1;

	  skip ();
	  if (ch != '(')
	    new_error ("mpfr_factor 1");

	  skip ();
	  mpfr_expr (tmp);

	  if (ch != ')')
	    new_error ("mpfr_factor 2");

	  skip ();

	  (int_fns[i].evalfn) (ret, tmp);

          mpfr_clear (tmp);

	  return;
	}
    }

  for (int i = 0; twoarg_fns[i].spelling != 0; i++)
    {
      char *spelling = twoarg_fns[i].spelling;
      int len = strlen (spelling);
      if (strncmp (spelling, bp - 1, len) == 0 && ! isalnum (bp[-1 + len]))
	{
	  bp += len - 1;

	  skip ();
	  if (ch != '(')
	    new_error ("mpfr_factor 1");

	  skip ();
	  mpfr_expr (tmp);

	  skip ();
	  if (ch != ',')
	    new_error ("mpfr_factor 2");

	  skip ();
	  mpfr_expr (tmp2);

	  skip ();
	  if (ch != ')')
	    new_error ("mpfr_factor 3");

	  skip ();

	  (twoarg_fns[i].evalfn) (ret, tmp, tmp2, default_rmode);

          mpfr_clear (tmp);
          mpfr_clear (tmp2);

	  return;
	}
    }

  for (int i = 0; int_arg_fns[i].spelling != 0; i++)
    {
      char *spelling = int_arg_fns[i].spelling;
      int len = strlen (spelling);
      if (strncmp (spelling, bp - 1, len) == 0 && ! isalnum (bp[-1 + len]))
	{
	  long int n;

	  bp += len - 1;

	  skip ();
	  if (ch != '(')
	    new_error ("mpfr_factor 1");

	  skip ();
	  mpfr_expr (tmp);
	  n = mpfr_get_ui (tmp, default_rmode);

	  skip ();
	  if (ch != ',')
	    new_error ("mpfr_factor 2");

	  skip ();
	  mpfr_expr (tmp2);

	  skip ();
	  if (ch != ')')
	    new_error ("mpfr_factor 3");

	  skip ();

	  (int_arg_fns[i].evalfn) (ret, n, tmp2, default_rmode);

          mpfr_clear (tmp);
          mpfr_clear (tmp2);

	  return;
	}
    }

  for (int i = 0; log_fns[i].spelling != 0; i++)
    {
      char *spelling = log_fns[i].spelling;
      int len = strlen (spelling);
      if (strncmp (spelling, bp - 1, len) == 0 && ! isalnum (bp[-1 + len]))
	{
	  int sign;

	  bp += len - 1;

	  skip ();
	  if (ch != '(')
	    new_error ("mpfr_factor 1");

	  skip ();
	  mpfr_expr (tmp);

	  skip ();
	  if (ch != ')')
	    new_error ("mpfr_factor 2");

	  skip ();

	  (log_fns[i].evalfn) (ret, &sign, tmp, default_rmode);

          mpfr_clear (tmp);

	  return;
	}
    }

  if (ch == '(')
    {
      skip ();
      mpfr_expr (ret);
      if (ch == ')')
	skip ();
      else
	new_error ("mpfr_factor 3");
    }
  else
    mpfr_number (ret);
  if (ch == '!')
    {
      if (mpfr_integer_p (ret) == 0 || mpfr_fits_ulong_p (ret, GMP_RNDN) == 0)
	new_error ("mpfr_factor 4");
      mpfr_fac_ui (ret, mpfr_get_ui (ret, GMP_RNDN), default_rmode);
      skip ();
    }

   mpfr_clear (tmp);
}

static void
mpfr_number (mpfr_t ret)
{
  char *endp;

  /* Exit */
  if (strncmp ("exit", bp - 1, 4) == 0 || strncmp ("x", bp - 1, 1) == 0
	|| strncmp ("quit", bp - 1, 4) == 0 || strncmp ("q", bp - 1, 1) == 0)
    exit(0);

  /* Pi */
  if (strncmp ("pi", bp - 1, 2) == 0 && ! isalnum (bp[1]))
    {
      bp += 2 - 1;
      skip ();
      mpfr_const_pi (ret, default_rmode);
      return;
    }
  /* Euler's constant */
  if (strncmp ("eu", bp - 1, 2) == 0 && ! isalnum (bp[1]))
    {
      bp += 2 - 1;
      skip ();
      mpfr_const_euler (ret, default_rmode);
      return;
    }
  // Catalan's constant
  if (strncmp ("ca", bp - 1, 2) == 0 && ! isalnum (bp[1]))
    {
      bp += 2 - 1;
      skip ();
      mpfr_const_catalan (ret, default_rmode);
      return;
    }
  if (ch == '$') /* previous result */
    {
      if (! previous_result_valid_flag)
	new_error ("mpfr_number 1");
      skip ();
      mpfr_set (ret, previous_result, default_rmode);
      return;
    }
  if (ch != '.' && (ch < '0' || ch > '9'))
    new_error ("mpfr_number 2");
  mpfr_strtofr (ret, bp - 1, &endp, 10, default_rmode);
  if (endp == bp - 1)
    new_error ("mpfr_number 3");
  bp = endp;
  skip ();
}

int
main (int argc, char *argv[])
{
  int nl_flag = 1;
  int prec = 128;

#define ZIV_PREC 32
  mpfr_t tmp;
  mpfr_t tmp_long[2]; /* tmp_long[0].prec = prec + ZIV_PREC */
                      /* tmp_long[1].prec = prec + 2 * ZIV_PREC */
  mp_exp_t tmp_long_exp[2];
//  char tmp_long_str[2][2500]; /* ceil(8192 * log10(2)) = 2467 */
  char *tmp_long_str[2] = {0, 0}; /* ceil(8192 * log10(2)) = 2467 */

  while (argc >= 2)
    {
      if (!strcmp (argv[1], "-n"))
	{
	  nl_flag = 0;
	  argv++;
	  argc--;
	}
      else if (!strcmp (argv[1], "-prec"))
	{
	  prec = atoi (argv[2]);

          /* fix!: 2011-10-05 by T.Kouya */
	  tmp_long_str[0] = (char *)malloc(sizeof(char) * ceil((prec + ZIV_PREC) * log10(2)));
	  tmp_long_str[1] = (char *)malloc(sizeof(char) * ceil((prec + ZIV_PREC) * log10(2)));
          if((prec > 1048576) || (tmp_long_str[0] == NULL) || (tmp_long_str[1] == NULL)) { printf("Precision is too large!(%d)\n", prec); exit(0);} 

	  argv += 2;
	  argc -= 2;
	}
      else if (!strcmp (argv[1], "-help") || !strcmp (argv[1], "-h"))
	{
          printf ("usage: %s [options] expr [expr ... ]\n", argv[0]);
	  printf ("   options:    -n      -- suppress newline\n");
	  printf ("               -prec n -- print n digits\n");
	  printf ("               -help   -- make demons drop from your nose\n");
	  exit (0);
	}
      else
	break;
    }
  if (tmp_long_str[0] == 0)
    tmp_long_str[0] = (char *)malloc(sizeof(char) * ceil((prec + ZIV_PREC) * log10(2)));
  if (tmp_long_str[1] == 0)
    tmp_long_str[1] = (char *)malloc(sizeof(char) * ceil((prec + ZIV_PREC) * log10(2)));

  if (argc >= 2)
    {
      int i;

      mpfr_init2 (tmp, prec);
      mpfr_init2 (tmp_long[0], prec + ZIV_PREC);
      mpfr_init2 (tmp_long[1], prec + ZIV_PREC * 2);
      mpfr_init2 (previous_result, prec);

      for (i = 1; i < argc; i++)
	{
	  expr_start_p = argv[i];
	  expr_end_p = expr_end_p + strlen (expr_start_p);
	  bp = expr_start_p;
	  skip ();
	  if (setjmp (top) == 0)
	    {
	      mpfr_expr (tmp_long[0]); /* only one exec! */
	      mpfr_get_str (tmp_long_str[0], &tmp_long_exp[0], 10,
                            ceil(prec * log10(2.0)), tmp_long[0], default_rmode);
	      mpfr_get_str (tmp_long_str[1], &tmp_long_exp[1], 10,
                            ceil(prec * log10(2.0)), tmp_long[0], GMP_RNDZ);
              mpfr_set (tmp, tmp_long[0], default_rmode);
	      mpfr_out_str (stdout, 10, 0, tmp, default_rmode);
	      if (nl_flag)
		puts ("");
	      mpfr_set (previous_result, tmp, default_rmode);
	      previous_result_valid_flag = 1;
	    }
	}
    }
  else
    {
#define BUFSIZE 1024
      char buf[BUFSIZE];

      mpfr_init2 (tmp, prec);
      mpfr_init2 (tmp_long[0], prec + ZIV_PREC);
      mpfr_init2 (tmp_long[1], prec + ZIV_PREC * 2);
      mpfr_init2 (previous_result, prec);

      for (;;)
	{
	  fputs ("eval> ", stdout);
	  bp = fgets (buf, BUFSIZE, stdin);
	  if (strncmp (bp, "\n", 1) == 0)
	    break;
	  skip ();
	  if (setjmp (top) == 0)
	    {
	      mpfr_expr (tmp_long[0]); /* only one exec! */
	      mpfr_get_str (tmp_long_str[0], &tmp_long_exp[0], 10,
                            ceil(prec * log10(2.0)), tmp_long[0], default_rmode);
	      mpfr_get_str (tmp_long_str[1], &tmp_long_exp[1], 10,
                            ceil(prec * log10(2.0)), tmp_long[0], GMP_RNDZ);
              mpfr_set (tmp, tmp_long[0], default_rmode);
              mpfr_out_str (stdout, 10, 0, tmp, default_rmode);
	      if (nl_flag)
		puts ("");
	      mpfr_set (previous_result, tmp, default_rmode);
	      previous_result_valid_flag = 1;
	    }
	}
    }

  mpfr_clear (tmp);
  mpfr_clear (tmp_long[0]);
  mpfr_clear (tmp_long[1]);
  mpfr_clear (previous_result);

  exit (0);
}


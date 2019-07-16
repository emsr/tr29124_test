{demo program from mp_intro, with Delphi use -cc command line switch}
program mp_intro;

uses
  mp_types, mp_base, mp_numth;
var
  a: mp_int;
begin
  writeln('MPArith V', MP_VERSION, ' with MAXDigits = ',MAXDigits, ', MP_MAXBIT = ',MP_MAXBIT);
  mp_init(a);
  mp_fib(271,a);
  writeln('Fibonacci(271) = ',mp_decimal(a));
  mp_2expt(a,MP_MAXBIT-1);
  writeln('Number of decimal digits of 2^',MP_MAXBIT-1,' = ', mp_radix_size(a,10));
  mp_clear(a);
end.

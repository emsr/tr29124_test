{Main test program for MPArith/mp_int,  (c) W.Ehrhardt 2004-2012}
{Many test cases derived from M.J. Fromberger's IMATH           }
{Other tests are generated with ARIBAS, Pari/GP, NX and others  }

program t_check;

{$x+}  {pchar I/O}

{$i STD.INC}
{$i mp_conf.inc}

{$i+}  {RTE on I/O error}

{$ifndef FPC}
{$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  BTypes, mp_types, mp_base,
  {$ifdef MPC_Diagnostic}
    mp_supp,
  {$endif}
  mp_prime, mp_modul, mp_numth, mp_pfu, mp_prng, mp_rsa, pfdu, t_check1;


{---------------------------------------------------------------------------}
procedure gcd32_test;
type
  t3l= array[1..3] of longint;
const
  NTC = 20;
  testcase: array[1..NTC] of t3l =
             (  (0, 0, 0),
                (0, 1, 1),
                (0, 2, 2),
                (0, 3, 3),
                (0, 4, 4),
                (1, 1, 1),
                (1, 2, 1),
                (1, 3, 1),
                (1, 4, 1),
                (2, 2, 2),
                (2, 3, 1),
                (2, 4, 2),
                (3, 3, 3),
                (3, 4, 1),
                (6,12, 6),
                (444444444, 666666666, 222222222),
                ($40000000, $60000000, $20000000),
                (longint($80000000), longint($C0000000), $40000000),
                (longint($C6AEA155), longint($84746B8E), 1111111111),
                (MaxLongint, -1, 1));
var
  i,cnt,f: integer;
  TA, TB, TC: longint;
begin
  writeln('Test gcd32');
  cnt := 0;
  f := 0;
  for i:=1 to NTC do begin
    inc(cnt);
    TA := testcase[i][1];
    TB := testcase[i][2];
    TC := testcase[i][3];
    if gcd32u(TA, TB) <> TC then begin
      writeln('gcd32u: error for testcase: ',i);
      inc(f);
    end;
    inc(cnt);
    if gcd32u(TB, TA) <> TC then begin
      writeln('gcd32u: error for swapped testcase: ',i);
      inc(f);
    end;
    if (TA>=0) and (TB>=0) then begin
      inc(cnt);
      if gcd32(TA, TB) <> TC then begin
        writeln('gcd32: error for (++) testcase : ',i);
        inc(f);
      end;
      inc(cnt);
      if gcd32(TA, -TB) <> TC then begin
        writeln('gcd32: error for (+-) testcase : ',i);
        inc(f);
      end;
      inc(cnt);
      if gcd32(-TA, TB) <> TC then begin
        writeln('gcd32: error for (-+) testcase : ',i);
        inc(f);
      end;
      inc(cnt);
      if gcd32(-TA, -TB) <> TC then begin
        writeln('gcd32: error for (--) testcase : ',i);
        inc(f);
      end;
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure xgcd32_test;
var
  a,b,u,v,d,c,i: longint;
{$ifdef HAS_INT64}
  a64,b64: int64;
{$else}
  a64,b64: comp;
{$endif}
  f: integer;
begin
  writeln('Test xgcd32');
  f := 0;
  for a:=-100 to 100 do begin
    if f<>0 then break;
    for b:=-100 to 100 do begin
      if f<>0 then break;
      xgcd32(a,b,u,v,d);
      c := gcd32(a,b);
      if c<>d then begin
        writeln('GCD diff: a=',a,',  b=',b,',  d=',d, ',  c=',c);
        inc(f);
      end;
      {Here there will be no overflows with |a|,|b| <= 100}
      if a*u+b*v<>d  then begin
        writeln('a*u+b*v<>d for a=',a,',  b=',b,',  u=',v, ',  v=',v, ',  d=',d);
        inc(f);
      end;
    end;
  end;

  for i:=1 to 50000 do begin
    if i and $3FF = 0 then write('.');
    a := mp_random_int;
    b := mp_random_int;
    xgcd32(a,b,u,v,d);
    c := gcd32(a,b);
    if c<>d then begin
      writeln('GCD diff: a=',a,',  b=',b,',  d=',d, ',  c=',c);
      inc(f);
      break;
    end;
    a64 := a;
    b64 := b;
    if a64*u + b64*v <> d  then begin
      writeln('a*u+b*v<>d for a=',a,',  b=',b,',  u=',v, ',  v=',v, ',  d=',d);
      inc(f);
      break;
    end;
  end;
  writeln;
  writeln('No. of checks failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure IsLong_Test;
var
  i,cnt,k,f: word;
  xl: longint;
  isli, is32: boolean;
begin
  f := 0;
  cnt := 0;
  mp_zero(a);
  writeln('Test mp_is_longint: ');
  for i:=0 to 64 do begin
    isli := mp_is_longint(a,xl);
    is32 := mp_bitsize(a)<32;
    if isli<>is32 then inc(f)
    else begin
      if isli then begin
        mp_set_int(c,xl);
        if mp_is_ne(a,c) then inc(f);
      end;
    end;
    mp_chs(a,b);
    isli := mp_is_longint(b,xl);
    is32 := mp_bitsize(b)<32;
    if isli<>is32 then inc(f)
    else begin
      if isli then begin
        mp_set_int(d,xl);
        if mp_is_ne(b,d) then inc(f);
      end;
    end;
    inc(cnt,2);
    mp_shl(a,1,a);
    mp_inc(a);
  end;
  for i:=0 to 50 do begin
    write('.');
    for k:=1 to 1000 do begin
      mp_rand_ex(a, 1 + i div DIGIT_BIT, false);
      if random(2)=0 then mp_chs(a,a);
      isli := mp_is_longint(a,xl);
      is32 := mp_bitsize(a)<32;
      if isli<>is32 then inc(f)
      else begin
        if isli then begin
          mp_set_int(c,xl);
          if mp_is_ne(a,c) then inc(f);
        end;
      end;
      inc(cnt);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure TestPrimes16Index;
var
  n: word;
  i: integer;
  OK: boolean;
begin
  write('Test Primes16Index ... ');
  OK := true;
  for i:=1 to NumPrimes16 do begin
    if i <> Primes16Index(Primes16[i]) then begin
      OK := false;
      inc(totalfailed);
      writeln('Unstable for i=',i);
      continue;
    end;
  end;
  for n:=0 to $FFFF do begin
    i := Primes16Index(n);
    if i>NumPrimes16 then begin
      OK := false;
      inc(totalfailed);
      writeln('Index too large, n=',n, ',  i=',i);
      continue;
    end;
    if i=0 then begin
      OK := false;
      inc(totalfailed);
      writeln('Index 0, n=',n);
      continue;
    end;
    if (n>1) and (Primes16[i]>n) then begin
      OK := false;
      inc(totalfailed);
      writeln('Primes16[i] > n, n=',n);
    end;
    if (i<NumPrimes16) and (Primes16[i+1]<=n) then begin
      OK := false;
      inc(totalfailed);
      writeln('Primes16[i] not maximal, n=',n, ',  i=',i);
    end;
  end;
  if OK then writeln('passed');
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_pi_lphi;
type
  ttpprec = record x,a,phi,ppi: longint end;
const
  NTC = 13;
  tcases: array[1..NTC] of ttpprec = (
            (x:     10000; a:  8;  phi:      1711;  ppi:     1229), {Riesel [40]}
            (x:    100000; a:  7;  phi:     18053;  ppi:     9592), {...}
            (x:  15485864; a: 53;  phi:   1568715;  ppi:  1000000), {Lehmer [39]}
            (x:  20000000; a: 58;  phi:   1988057;  ppi:  1270607), {...}
            (x:  25000000; a: 61;  phi:   2458751;  ppi:  1565927),
            (x:  32452845; a: 66;  phi:   3140783;  ppi:  2000000),
            (x:  33000000; a: 66;  phi:   3193726;  ppi:  2031667),
            (x:  37000000; a: 67;  phi:   3569832;  ppi:  2261623),
            (x:  40000000; a: 68;  phi:   3847872;  ppi:  2433654),
            (x:  90000000; a: 25;  phi:  10826326;  ppi:  5216954),
            (x: 100000000; a: 90;  phi:   9110819;  ppi:  5761455),
            (x: 999000000; a: 40;  phi: 107108994;  ppi: 50799577),
            (x:1000000000; a: 40;  phi: 107216231;  ppi: 50847534));
var
  i,cnt,f: integer;
  tphi, tpi: longint;
begin
  write('Test primepi32 / lsumphi32 ');
  cnt := 0;
  f := 0;
  for i:=1 to NTC do with tcases[i] do begin
    write('.');
    tpi  := primepi32(x);
    tphi := lsumphi32(x,a);
    inc(cnt,2);
    if tpi<>ppi then begin
      inc(f);
      writeln('Failure primepi32 for i =', i);
    end;
    if tphi<>phi then begin
      inc(f);
      writeln('Failure lphi for i =', i);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_prime32;
const  {http://primes.utm.edu/nthprime/index.php#nth}
  p10: array[0..8] of longint = (2,29,541,7919,104729,1299709,15485863,
                                 179424673,2038074743);
var
  i,cnt,f: integer;
  p,k: longint;
begin
  write('Test prime32 ... ');
  cnt := 0;
  f   := 0;
  k   := 1;
  for i:=0 to 8 do begin
    if i<>0 then k := 10*k;
    p := prime32(k);
    inc(cnt);
    if p<>p10[i] then begin
      inc(f);
      writeln('Failure for k =', k);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_4sq;
{$ifdef BIT16}
const
  BITSI = 4;
{$else}
const
  BITSI = 8;
{$endif}
var
  i,bits,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  write('mp_4sq test ');
  for i:=1 to 100 do begin
    if odd(i) then bits := 50*BITSI
    else begin
      bits := i*BITSI;
      write('.');
    end;
    mp_rand_bits_ex(e,bits,false);
    if random(4)=0 then mp_shl(e,random(10),e);
    case random(3) of
      0: mp_4sq(e,a,b,c,d);
      1: mp_4sq_sa(e,a,b,c,d);
      2: mp_4sq(e,a,b,c,d);
    end;
    inc(cnt);
    mp_sqr(a,s);
    mp_sqr(b,t); mp_add(s,t,s);
    mp_sqr(c,t); mp_add(s,t,s);
    mp_sqr(d,t); mp_add(s,t,s);
    if mp_is_ne(s,e) then begin
      mp_writeln('e = ',e);
      mp_writeln('s = ',s);
      mp_writeln('a = ',a);
      mp_writeln('b = ',b);
      mp_writeln('c = ',c);
      mp_writeln('d = ',d);
      inc(f);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{$ifndef BIT16}
{---------------------------------------------------------------------------}
procedure astr_test;
var
  i,cnt,f: integer;
  bits: longint;
  r: word;
  s: ansistring;
begin
  f := 0;
  cnt := 0;
  writeln('mp_(read)_radix_astr test: ');
  for i:=1 to 100 do begin
    if i and 2 = 0 then write('.');
    r := 2+random(MaxRadix-1);
    mp_show_plus := random(2)=0;
    bits := sqr(longint(i))*10;
    mp_rand_bits_ex(a,bits,false);
    if random(2)=0 then mp_chs(a,a);
    s := mp_radix_astr(a,r);
    mp_read_radix_astr(b,s,r);
    inc(cnt);
    if mp_is_ne(a,b) then begin
      writeln('Failed for i/bits/radix=',i,'/',bits,'/',r);
      inc(f);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;
{$endif}


{---------------------------------------------------------------------------}
procedure test_catalan;
var
  i,cnt,f: integer;
  p: word;
  n,am,bm: longint;

{$ifdef MP_16BIT}
  {$ifdef BIT16}
    const imax=30;
  {$else}
    const imax=35;
  {$endif}
{$else}
  const imax=40;
{$endif}
const
  NSC = 19;
  sc: array[0..NSC] of longint =  {http://oeis.org/A000108}
        (1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012,
         742900, 2674440, 9694845, 35357670, 129644790, 477638700, 1767263190);
begin
  f := 0;
  cnt := 0;
  writeln('mp_catalan test: ');
  for i:=0 to NSC do begin
    inc(cnt);
    mp_catalan(i,a);
    mp_sub_int(a,sc[i],a);
    if not mp_is0(a) then begin
      inc(f);
      writeln('Small check failed for i=',i);
    end;
  end;

  for i:=2 to imax do begin
    inc(cnt);
    {Check: (1/2)*Catalan((p^2-1)/2) mod p^2 = 2^p-1 mod p^2}
    p := primes16[2*i];
    n := sqr(longint(p));
    mp_catalan((n-1) div 2, a);
    mp_shr1(a);
    mp_mersenne(p,b);
    mp_mod_int(a,n,am);
    mp_mod_int(b,n,bm);
    if am<>bm then begin
      inc(f);
      writeln('Failed for p=',p);
    end;
    if odd(i) then write('.');
  end;

  {check C(n) = sum(C(i)*C(n-1-i), i=0..n-1}
  mp_catalan(1000,a);
  for i:=0 to 999 do begin
    if i and 63 = 0 then write('+');
    mp_catalan(i,b);
    mp_catalan(999-i,c);
    mp_mul(b,c,d);
    mp_sub(a,d,a);
  end;
  inc(cnt);
  if not mp_is0(a) then begin
    inc(f);
    writeln('Summation check failed');
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_poch_perm;
var
  i,n,cnt,f: integer;
const
  tpoch : array[-6..6, 0..7] of longint = (
            ( 1,  -6,  30,  -120,   360,  -720,    720,       0 ),
            ( 1,  -5,  20,   -60,   120,  -120,      0,       0 ),
            ( 1,  -4,  12,   -24,    24,     0,      0,       0 ),
            ( 1,  -3,   6,    -6,     0,     0,      0,       0 ),
            ( 1,  -2,   2,     0,     0,     0,      0,       0 ),
            ( 1,  -1,   0,     0,     0,     0,      0,       0 ),
            ( 1,   0,   0,     0,     0,     0,      0,       0 ),
            ( 1,   1,   2,     6,    24,   120,    720,    5040 ),
            ( 1,   2,   6,    24,   120,   720,   5040,   40320 ),
            ( 1,   3,  12,    60,   360,  2520,  20160,  181440 ),
            ( 1,   4,  20,   120,   840,  6720,  60480,  604800 ),
            ( 1,   5,  30,   210,  1680, 15120, 151200, 1663200 ),
            ( 1,   6,  42,   336,  3024, 30240, 332640, 3991680 ));
const
  tperm : array[0..8, 0..8] of word = (
            ( 1, 0,  0,   0,    0,    0,     0,     0,     0 ),
            ( 1, 1,  0,   0,    0,    0,     0,     0,     0 ),
            ( 1, 2,  2,   0,    0,    0,     0,     0,     0 ),
            ( 1, 3,  6,   6,    0,    0,     0,     0,     0 ),
            ( 1, 4, 12,  24,   24,    0,     0,     0,     0 ),
            ( 1, 5, 20,  60,  120,  120,     0,     0,     0 ),
            ( 1, 6, 30, 120,  360,  720,   720,     0,     0 ),
            ( 1, 7, 42, 210,  840, 2520,  5040,  5040,     0 ),
            ( 1, 8, 56, 336, 1680, 6720, 20160, 40320, 40320 ));
{$ifdef BIT16}
const
  ifmax = 8;
{$else}
const
  ifmax = 14;
{$endif}

begin
  f := 0;
  cnt := 0;
  write('Test mp_product / mp_poch / mp_perm ');

  for n:=-6 to 6 do begin
    write('.');
    for i:=0 to 7 do begin
      inc(cnt);
      mp_poch(n,i,a);
      mp_sub_int(a,tpoch[n,i],b);
      if not mp_iszero(b) then begin
        writeln('mp_poch tab fail n,i =',n,',',i);
        inc(f);
      end;
    end;
  end;
  inc(cnt,2);
  mp_read_decimal_str(b,'6729457515998774705268610104629483525487942150402731785206632546304000000000');
  mp_poch(60,40,a);
  if mp_is_ne(a,b) then begin
    writeln('mp_poch(60,40) failed');
    inc(f);
  end;
  mp_read_decimal_str(b,'45159067006671175434518243598614826795313039459718126370816000000000');
  mp_poch(-70,40,a);
  if mp_is_ne(a,b) then begin
    writeln('mp_poch(-70,40) failed');
    inc(f);
  end;

  for n:=0 to 8 do begin
    write('.');
    for i:=0 to 8 do begin
      inc(cnt);
      mp_perm(n,i,a);
      mp_sub_int(a,tperm[n,i],b);
      if not mp_iszero(b) then begin
        writeln('mp_perm tab fail n,i = ',n,',',i);
        inc(f);
      end;
    end;
  end;
  inc(cnt);
  mp_read_decimal_str(b,'3420189997285434351598365756849270130350567547616624640000000000');
  mp_perm(60,40,a);
  if mp_is_ne(a,b) then begin
    writeln('mp_perm(60,40) failed');
    inc(f);
  end;

  for i:=0 to ifmax do begin
    write('.');
    inc(cnt);
    n := 1000*i;
    mp_fact(n,a);
    mp_product(b,1,n);
    if mp_is_ne(a,b) then begin
      writeln('mp_product failed for i = ',i);
      inc(f);
    end;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_EulerPhi;
const
  NE1=200;
  ea: array[0..NE1] of byte = {with(numtheory): [seq(phi(i),i=0..200)];}
        (0,1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8,
        12,10,22,8,20,12,18,12,28,8,30,16,20,16,24,12,36,
        18,24,16,40,12,42,20,24,22,46,16,42,20,32,24,52,
        18,40,24,36,28,58,16,60,30,36,32,48,20,66,32,44,
        24,70,24,72,36,40,36,60,24,78,32,54,40,82,24,64,
        42,56,40,88,24,72,44,60,46,72,32,96,42,60,40,100,
        32,102,48,48,52,106,36,108,40,72,48,112,36,88,56,
        72,58,96,32,110,60,80,60,100,36,126,64,84,48,130,
        40,108,66,72,64,136,44,138,48,92,70,120,48,112,72,
        84,72,148,40,150,72,96,60,120,48,156,78,104,64,132,
        54,162,80,80,82,166,48,156,64,108,84,172,56,120,80,
        116,88,178,48,180,72,120,88,144,60,160,92,108,72,
        190,64,192,96,96,84,196,60,198,80);
const
  NE2=10;
  nl: array[1..NE2] of longint = ({random values and 2^30}
        35690018,1809741004,112970497,1577004683,1997451014,
        1419616041,862708834,403211554,2113914896,$40000000);
  el: array[1..NE2] of longint = (
        13374720,904870500,112970496,1508439240,920131488,
        946148112,398173152,201605776,1056571776,$20000000);
var
  i,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  writeln('Test EulerPhi32 ...');
  for i:=0 to NE1 do begin
    inc(cnt);
    if EulerPhi32(i) <> ea[i] then begin
      writeln('Failed for i = ',i);
      inc(f);
    end;
  end;
  for i:=1 to NE2 do begin
    inc(cnt);
    if EulerPhi32(nl[i]) <> el[i] then begin
      writeln('Failed for ', nl[i]);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_carmichael;
const
  NC1=200;
  ca: array[0..NC1] of byte = {with(numtheory): [seq(lambda(n), n=0..200)];}
        (0,1,1,2,2,4,2,6,2,6,4,10,2,12,6,4,4,16,6,18,4,6,
         10,22,2,20,12,18,6,28,4,30,8,10,16,12,6,36,18,12,
         4,40,6,42,10,12,22,46,4,42,20,16,12,52,18,20,6,
         18,28,58,4,60,30,6,16,12,10,66,16,22,12,70,6,72,
         36,20,18,30,12,78,4,54,40,82,6,16,42,28,10,88,12,
         12,22,30,46,36,8,96,42,30,20,100,16,102,12,12,52,
         106,18,108,20,36,12,112,18,44,28,12,58,48,4,110,60,
         40,30,100,6,126,32,42,12,130,10,18,66,36,16,136,22,
         138,12,46,70,60,12,28,72,42,36,148,20,150,18,48,30,
         60,12,156,78,52,8,66,54,162,40,20,82,166,6,156,16,
         18,42,172,28,60,20,58,88,178,12,180,12,60,22,36,30,
         80,46,18,36,190,16,192,96,12,42,196,30,198,20);
const
  NC2=10;
  nl: array[1..NC2] of longint = ( {random values and 2^30}
        35690018,1809741004,112970497,1577004683,1997451014,
        1419616041,862708834,403211554,2113914896, $40000000);
  cl: array[1..NC2] of longint = (
        30960,452435250,112970496,68565420,12779604,
        236537028,99543288,201605776,33017868, $10000000);
var
  i,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  writeln('Test Carmichael32 ...');
  for i:=0 to NC1 do begin
    inc(cnt);
    if Carmichael32(i) <> ca[i] then begin
      writeln('Failed for i = ',i);
      inc(f);
    end;
  end;
  for i:=1 to NC2 do begin
    inc(cnt);
    if Carmichael32(nl[i]) <> cl[i] then begin
      writeln('Failed for ', nl[i]);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_is_carmichael;
var
  cnt,f: integer;
begin
  f := 0;
  cnt := 8;
  writeln('Test mp_is_pcarmichael ...');

  mp_read_decimal(a,'2140699681');
  if not mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_read_decimal(a,'2135039017');
  if mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_read_decimal(a,'725849475559');
  if mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_read_decimal(a,'9487582789999');
  if mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;
  mp_read_decimal(a,'4954039956700380001');
  if not mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_read_decimal(a,'3887636054124102392503405910694993617809');
  if not mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_read_decimal(a,'388764511027199999999999999414113587204032640000000000'+
                    '2943194650080514478879999999507163017775414944305441');
  if not mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  mp_nextprime(a);
  if mp_is_pcarmichael(a) then begin
    writeln('Failed for ', mp_decimal(a));
    inc(f);
  end;

  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_moebius;
const
  NM1=200;
  ma: array[1..NM1] of shortint = {with(numtheory): [[seq(mobius(x),x=1..200)];}
        (1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,
         -1,0,1,1,-1,0,0,1,0,0,-1,-1,-1,0,1,1,1,0,-1,
         1,1,0,-1,-1,-1,0,0,1,-1,0,0,0,1,0,-1,0,1,0,
         1,1,-1,0,-1,1,0,0,1,-1,-1,0,1,-1,-1,0,-1,1,
         0,0,1,-1,-1,0,0,1,-1,0,1,1,1,0,-1,0,1,0,1,
         1,1,0,-1,0,0,0,-1,-1,-1,0,-1,1,-1,0,-1,-1,1,
         0,-1,-1,1,0,0,1,1,0,0,1,1,0,0,0,-1,0,1,-1,
         -1,0,1,1,0,0,-1,-1,-1,0,1,1,1,0,1,1,0,0,-1,
         0,-1,0,0,-1,1,0,-1,1,1,0,1,0,-1,0,-1,1,-1,0,
         0,-1,0,0,-1,-1,0,0,1,1,-1,0,-1,-1,1,0,1,-1,
         1,0,0,-1,-1,0,-1,1,-1,0,-1,0,-1,0);
const
  NM2=10;
  nl: array[1..NM2] of longint = ({random values and 2^30}
        35690018,1809741004,112970497,1577004683,1997451014,
        1419616041,862708834,403211554,2113914896, $40000000);
  ml: array[1..NM2] of shortint = (-1, 0, -1, 1, 1, -1, 0, 1, 0, 0);
var
  i,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  writeln('Test Moebius32 ...');
  for i:=1 to NM1 do begin
    inc(cnt);
    if Moebius32(i) <> ma[i] then begin
      writeln('Failed for i = ',i);
      inc(f);
    end;
  end;
  for i:=1 to NM2 do begin
    inc(cnt);
    if Moebius32(nl[i]) <> ml[i] then begin
      writeln('Failed for ', nl[i]);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_primroot;
const
  NR1=200;
  ra: array[2..NR1] of byte = ( {with(numtheory): [[seq(primroot(x),x=2..200)];}
        1,2,3,2,5,3,0,2,3,2,0,2,3,0,0,3,5,2,0,0,7,5,0,2,7,2,0,2,0,
        3,0,0,3,0,0,2,3,0,0,6,0,3,0,0,5,5,0,3,3,0,0,2,5,0,0,0,3,2,
        0,2,3,0,0,0,0,2,0,0,0,7,0,5,5,0,0,0,0,3,0,2,7,2,0,0,3,0,0,
        3,0,0,0,0,5,0,0,5,3,0,0,2,0,5,0,0,3,2,0,6,0,0,0,3,0,0,0,0,
        11,0,0,2,7,0,0,2,0,3,0,0,0,2,0,0,7,0,0,3,0,2,0,0,7,0,0,0,
        5,0,0,2,0,6,0,0,0,0,0,5,3,0,0,0,5,2,0,0,5,5,0,2,0,0,0,2,0,
        0,0,0,3,2,0,2,0,0,0,0,0,0,0,0,0,19,0,5,5,0,0,2,0,3,0);
const
  NM2=10;
  nl: array[1..NM2] of longint = ( {random values}
        112970497, 403211554, 2113914896, 1799212189, 31134751,
        607067782, 917676758, 2031081119, 1505450647, 29179441);
  ml: array[1..NM2] of byte = (5, 3, 0, 2, 6, 11, 15, 7, 3, 13);
var
  i,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  writeln('Test primroot32 ...');
  for i:=2 to NR1 do begin
    inc(cnt);
    if primroot32(i) <> ra[i] then begin
      writeln('Failed for i = ',i);
      inc(f);
    end;
  end;
  for i:=1 to NM2 do begin
    inc(cnt);
    if primroot32(nl[i]) <> ml[i] then begin
      writeln('Failed for ', nl[i]);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_order32;
const
  N1 = 100;
  M1 = 101;
  O1 : array[1..N1] of byte = ( {[seq(order(i,101),i=1..100)];}
        1,100,100,50,25,10,100,100,50,4,100,100,50,10,100,25,
        10,100,25,50,50,50,50,25,25,100,100,100,100,50,25,
        20,50,100,100,5,25,100,20,100,20,100,50,20,50,100,
        50,100,50,100,100,25,100,25,100,25,20,25,100,20,100,
        20,100,50,10,100,100,25,20,50,25,100,100,100,100,50,
        50,25,25,25,25,50,100,5,50,100,5,25,100,100,4,25,
        100,100,5,50,25,100,100,2);
  N2 = 97;
  M2 = 13*17;
  O2 : array[0..N2] of byte = ( {[seq(order(i,13*17),i=0..97)];}
         0,1,24,48,12,16,48,48,8,24,48,48,16,0,16,24,6,
         0,4,24,48,4,48,48,48,8,0,16,48,48,12,16,24,
         12,0,3,24,48,4,0,16,48,24,24,16,48,48,4,48,
         24,12,0,0,8,48,12,48,16,48,24,8,48,48,48,4,
         0,8,12,0,6,8,48,12,16,48,48,24,8,0,16,48,12,
         48,8,12,0,4,24,48,12,16,0,16,24,24,48,16,48);
var
  i,n: longint;
  cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  writeln('Test order32 ...');
  for i:=low(O1) to high(O1) do begin
    inc(cnt);
    n := order32(i,M1);
    if n <> O1[i] then begin
      inc(f);
      writeln('O1 failed: ',i);
    end;
  end;
  for i:=low(O2) to high(O2) do begin
    inc(cnt);
    n := order32(i,M2);
    if n <> O2[i] then begin
      inc(f);
      writeln('O2 failed: ',i);
      writeln(i);
    end;
  end;
  for i:=2 to 2000 do begin
    inc(cnt,2);
    n := sqr(i);
    if 2<>order32(i,i+1) then begin
      inc(f);
      writeln('Order 2 failed: ',i);
    end;
    if 4<>order32(i,n+1) then begin
      inc(f);
      writeln('Order 4 failed: ',i);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_sigmak;
var
  i,cnt,f: integer;
const
  sigma_0: array[1..15] of byte = ( {http://oeis.org/A000005}
             1,2,2,3,2,4,2,4,3,4,2,6,2,4,4);

  sigma_1: array[1..15] of byte = ( {http://oeis.org/A000203}
             1,3,4,7,6,12,8,15,13,18,12,28,14,24,24);

  sigma_2: array[1..15] of word = ( { http://oeis.org/A001157}
             1,5,10,21,26,50,50,85,91,130,122,210,170,250,260);

  sigma_3: array[1..15] of word = ( { http://oeis.org/A001158}
             1,9,28,73,126,252,344,585,757,1134,1332,2044,2198,3096,3528);

  sigma_4: array[1..15] of word = ( { http://oeis.org/A001159}
             1,17,82,273,626,1394,2402,4369,6643,10642,14642,22386,28562,40834,51332);

  sigma_5: array[1..15] of longint = ( { http://oeis.org/A001160}
             1,33,244,1057,3126,8052,16808,33825,59293,103158,161052,257908,371294,554664,762744);

const {Pari: sigma(n,k)}
  s11 = '2017807491561102999631484340618436673467747480392679929996';
  s12 = '2911221654831869466363023287520284782191148555791577871102530352168390412613776';
  s13 = '1547729197788406436578232909486595843792561058208796001316378245478351358596333115611712045580793104620576';
  s17 = '3595461272557299653174091079336363772054784552910504170610996694917'
       +'45365420254055103594069576816820962820097287171393649148064062176676896';
  s42 = '150130937545330708181211458454834608908702699731083721027349152152500';

begin
  f := 0;
  cnt := 0;
  writeln('Test mp_sigmak ...');
  for i:=1 to 15 do begin
    inc(cnt);
    mp_sigmak(0,i,a);
    if mp_cmp_int(a, sigma_0[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_0: ',i);
    end;
    inc(cnt);
    mp_sigmak(1,i,a);
    if mp_cmp_int(a, sigma_1[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_1: ',i);
    end;
    inc(cnt);
    mp_sigmak(2,i,a);
    if mp_cmp_int(a, sigma_2[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_2: ',i);
    end;
    inc(cnt);
    mp_sigmak(3,i,a);
    if mp_cmp_int(a, sigma_3[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_3: ',i);
    end;
    inc(cnt);
    mp_sigmak(4,i,a);
    if mp_cmp_int(a, sigma_4[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_4: ',i);
    end;
    inc(cnt);
    mp_sigmak(5,i,a);
    if mp_cmp_int(a, sigma_5[i])<>MP_EQ then begin
      inc(f);
      writeln('Fail sigma_5: ',i);
    end;
  end;

  mp_sigmak(11, 162000,a);
  mp_read_decimal_str(b,s11);
  inc(cnt);
  if mp_is_ne(a,b) then begin
    inc(f);
    writeln('S11 failed');
  end;

  mp_sigmak(12, 3456789,a);
  mp_read_decimal_str(b,s12);
  inc(cnt);
  if mp_is_ne(a,b) then begin
    inc(f);
    writeln('S12 failed');
  end;

  mp_sigmak(13, 123456789,a);
  mp_read_decimal_str(b,s13);
  inc(cnt);
  if mp_is_ne(a,b) then begin
    inc(f);
    writeln('S13 failed');
  end;

  mp_sigmak(17, 123456789,a);
  mp_read_decimal_str(b,s17);
  inc(cnt);
  if mp_is_ne(a,b) then begin
    inc(f);
    writeln('S17 failed');
  end;

  mp_sigmak(42, 42,a);
  mp_read_decimal_str(b,s42);
  inc(cnt);
  if mp_is_ne(a,b) then begin
    inc(f);
    writeln('S42 failed');
  end;

  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


(*
gp > polrootsmod(2*x^2 + 3*x + 0, 5)
%5 = [Mod(0, 5), Mod(1, 5)]~

gp > polrootsmod(2*x^2 + 3*x + 1, 5)
%6 = [Mod(2, 5), Mod(4, 5)]~

gp > polrootsmod(2*x^2 + 3*x + 2, 5)
%7 = []~

gp > polrootsmod(2*x^2 + 3*x + 3, 5)
%8 = [Mod(3, 5)]~

gp > polrootsmod(2*x^2 + 3*x + 4, 5)
%9 = []~

gp > polrootsmod(15*x^2 + 3*x + 1, 5)
%4 = [Mod(3, 5)]~

gp > p=2^89-1
%1 = 618970019642690137449562111

gp > polrootsmod(2*x^2 + 3*x + 77371252455336267181195265, p)
%2 = [Mod(154742504910672534362390527, 618970019642690137449562111)]~

gp > polrootsmod(2*x^2 + 3*x + 0, p)
%3 = [Mod(0, 618970019642690137449562111),
      Mod(309485009821345068724781054, 618970019642690137449562111)]~

gp > polrootsmod(2*x^2 + 3*x + 1, p)
%4 = [Mod(309485009821345068724781055, 618970019642690137449562111),
      Mod(618970019642690137449562110, 618970019642690137449562111)]~

gp > polrootsmod(2*x^2 + 3*x + 5, p)
%6 = []~

gp > polrootsmod(2*x^2 + 3*x + 11, p)
%8 = [Mod( 82970829808510104070914033, 618970019642690137449562111),
      Mod(226514180012834964653867021, 618970019642690137449562111)]~
*)

{---------------------------------------------------------------------------}
procedure test_squad_mod;
var
  i,k,cnt,f: integer;
  lf: boolean;
const
  qn: array[0..4] of byte = (2,2,0,1,0);
  x1: array[0..4] of byte = (0,2,0,3,0);
  x2: array[0..4] of byte = (1,4,0,0,0);
begin
  f := 0;
  cnt := 0;
  writeln('Test mp_squad_mod ...');
  mp_set(z,5);
  mp_set(a,2);
  mp_set(b,3);
  for k:=0 to 4 do begin
    inc(cnt);
    mp_set_int(c,k);
    i := mp_squad_mod(a,b,c,z,s,t);
    lf := false;
    if i=qn[k] then begin
      if i=1 then lf := mp_cmp_d(s,x1[k])<>MP_EQ
      else if i=2 then lf := (mp_cmp_d(s,x1[k])<>MP_EQ) or (mp_cmp_d(t,x2[k])<>MP_EQ);
    end
    else lf := true;
    if lf then begin
      inc(f);
      writeln('Failed for k=',k);
    end;
  end;

  mp_mersenne(89,z);

  inc(cnt);
  mp_read_decimal(c,'77371252455336267181195265');
  mp_read_decimal(d,'154742504910672534362390527');
  i := mp_squad_mod(a,b,c,z,s,t);
  if (i<>1) or mp_is_ne(s,d) then begin
    inc(f);
    writeln('Fail big 1');
  end;

  inc(cnt);
  mp_zero(c);
  mp_read_decimal(e,'309485009821345068724781054');
  i := mp_squad_mod(a,b,c,z,s,t);
  if (i<>2) or (not mp_is0(s)) or mp_is_ne(e,t) then begin
    inc(f);
    writeln('Fail big 2');
  end;

  inc(cnt);
  mp_set(c,1);
  mp_read_decimal(d,'309485009821345068724781055');
  mp_read_decimal(e,'618970019642690137449562110');
  i := mp_squad_mod(a,b,c,z,s,t);
  if (i<>2) or mp_is_ne(d,s) or mp_is_ne(e,t) then begin
    inc(f);
    writeln('Fail big 3');
  end;

  inc(cnt);
  mp_set(c,11);
  mp_read_decimal(d, '82970829808510104070914033');
  mp_read_decimal(e,'226514180012834964653867021');
  i := mp_squad_mod(a,b,c,z,s,t);
  if (i<>2) or mp_is_ne(d,s) or mp_is_ne(e,t) then begin
    inc(f);
    writeln('Fail big 4');
  end;

  inc(cnt);
  mp_set(c,5);
  i := mp_squad_mod(a,b,c,z,s,t);
  if i<>0 then begin
    inc(f);
    writeln('Fail big 5');
  end;

  {Test 2a not invertible}
  inc(cnt);
  mp_set(a,15);
  mp_set(b,3);
  mp_set(c,1);
  mp_set(z,5);
  i := mp_squad_mod(a,b,c,z,s,t);
  if (i<>1) or (mp_cmp_d(s,3)<>MP_EQ) then begin
    inc(f);
    writeln('Fail a=15');
  end;

  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_lucasuv;
var
  i,cnt,f: integer;
  j,n,nmax: longint;
const
  PN : array[0..20] of longint = ( {http://oeis.org/A000129} {Pell numbers}
         0, 1, 2, 5, 12, 29, 70, 169, 408, 985, 2378, 5741, 13860, 33461,
         80782, 195025, 470832, 1136689, 2744210, 6625109, 15994428);
const
  CP : array[0..20] of longint = ({http://oeis.org/A002203}  {Companion Pell numbers}
         2, 2, 6, 14, 34, 82, 198, 478, 1154, 2786, 6726, 16238, 39202, 94642,
         228486, 551614, 1331714, 3215042, 7761798, 18738638, 45239074);
const
  JU : array[0..20] of longint = ( {http://oeis.org/A001045} {Jacobsthal numbers}
         0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341, 683, 1365, 2731, 5461,
         10923, 21845, 43691, 87381, 174763, 349525);
const
  JL : array[0..20] of longint = ({http://oeis.org/A014551}  {Jacobsthal-Lucas numbers}
         2, 1, 5, 7, 17, 31, 65, 127, 257, 511, 1025, 2047, 4097, 8191,
         16385, 32767, 65537, 131071, 262145, 524287, 1048577);
const
  P128 = '3497379255757941172020851852070562919437964212608'; {Pari: a=[2,1;1,0]^128; a[1,2]}
  MAXMAX = 2000000;
begin
  f := 0;
  cnt := 0;
  write('Test mp_lucasuv ');
  nmax := MaxMersenne div 2;
  if nmax > MAXMAX then nmax := MAXMAX;

  {Pell numbers p=a=2, q=b=-1}
  mp_set_int(a,2);
  mp_set_int(b,-1);
  for i:=0 to 20 do begin
    if odd(i) then write('.');
    inc(cnt);
    mp_lucasuv(a,b,i,c,d);
    mp_lucasu(a,b,i,e);
    mp_lucasv(a,b,i,z);
    mp_set_int(s,PN[i]);
    mp_set_int(t,CP[i]);
    if mp_is_ne(c,s) or mp_is_ne(d,t) or mp_is_ne(e,s) or mp_is_ne(z,z) then begin
      inc(f);
      writeln('Pell numbers failed for i=',i);
    end;
  end;

  inc(cnt);
  mp_lucasu(a,b,128,e);
  mp_read_decimal(s,P128);
  if mp_is_ne(e,s) then begin
    inc(f);
    writeln('Pell(128) failed');
  end;

  {Jacobsthal numbers p=a=1, q=b=-2}
  mp_set_int(a,1);
  mp_set_int(b,-2);
  for i:=0 to 20 do begin
    if odd(i) then write('.');
    inc(cnt);
    mp_lucasuv(a,b,i,c,d);
    mp_lucasu(a,b,i,e);
    mp_lucasv(a,b,i,z);
    mp_set_int(s,JU[i]);
    mp_set_int(t,JL[i]);
    if mp_is_ne(c,s) or mp_is_ne(d,t) or mp_is_ne(e,s) or mp_is_ne(z,z) then begin
      inc(f);
      writeln('Jacobsthal failed for i=',i);
    end;
  end;
  n := 21;
  while n<=nmax do begin
    write('.');
    inc(cnt);
    mp_lucasuv(a,b,n,c,d);
    mp_2expt(s, n);
    if odd(n) then mp_inc(s) else mp_dec(s);
    mp_div_int(s,3,@s,j);
    mp_2expt(t, n);
    if odd(n) then mp_dec(t) else mp_inc(t);
    if mp_is_ne(c,s) or mp_is_ne(d,t) then begin
      inc(f);
      writeln('Jacobsthal failed for n=',n);
    end;
    n := 3*n+1;
  end;

  mp_set(a,1);
  mp_set_int(b,-1);
  n := 0;
  while n<=nmax do begin
    write('.');
    inc(cnt);
    mp_lucasuv(a,b,n,c,d);
    mp_fib(n,s);
    mp_lucas(n,t);
    if mp_is_ne(c,s) or mp_is_ne(d,t) then begin
      inc(f);
      writeln('Fib/Lucus failed for n=',n);
    end;
    n := 3*n+1;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_qnr;
var
  i,cnt,f: integer;
  j: longint;
type
 trec = record
          qnr: longint;
          ns : pchar8;
        end;
{Calculated with: CALC - a number theory calculator, K.R.Matthews}
{available from http://www.numbertheory.org/calc/krm_calc.html}
{and with the Pari script for composite and prime values:}
(*
qnr(n)={local(r); r=-1; if (issquare(n),r=-1, for(i=2, n-1, if (kronecker(r=i,n)==-1, break, r=-1))); r}
*)
const
  NT = 23;
  tests: array[1..NT] of trec = (
           (qnr: 7   ; ns: '1031'),   {prime}
           (qnr: 11  ; ns: '311'),
           (qnr: 47  ; ns: '3818929'),
           (qnr: 67  ; ns: '48473881'),
           (qnr: 73  ; ns: '120293879'),
           (qnr: 83  ; ns: '131486759'),
           (qnr: 17  ; ns: '9741602165436667319'),
           (qnr: 13  ; ns: '261020784828312504681720094716692999449'),
           (qnr: 67  ; ns: '281857014562698042780079535303369444209'),

           (qnr: 5   ; ns: '1023'),  {odd composite}
           (qnr: 7   ; ns: '15'),
           (qnr: 17  ; ns: '231'),
           (qnr: 41  ; ns: '25311'),
           (qnr: 67  ; ns: '268719'),
           (qnr: 103 ; ns: '472585905'),
           (qnr: 109 ; ns: '715236599'),
           (qnr: 61  ; ns: '925791654865'),
           (qnr: 19  ; ns: '246257702670521861564047901568556771129'),

           (qnr: 11  ; ns: '886565939000'), {even}
           (qnr: 13  ; ns: '602288042610'),
           (qnr: 17  ; ns: '1000992936060'),
           (qnr: 139 ; ns: '760636202914'),
           (qnr: 103 ; ns: '16506048235600416393176156735904')
        );
const
  q100: array[1..100] of shortint = ( {vector(100, n, qnr(n))}
          -1, -1, 2, -1, 2, -1, 3, 3, -1, 7, 2, 5, 2, 11, 7, -1, 3, 5, 2, 3, 2, 3,
           5, 13, -1, 3, 2, 3, 2, 7, 3, 3, 5, 7, 2, -1, 2, 5, 7, 7, 3, 5, 2, 7, 2, 3, 5, 5,
          -1, 3, 2, 5, 2, 13, 3, 11, 5, 5, 2, 7, 2, 5, 5, -1, 3, 7, 2, 3, 2, 3, 7, 5, 5,
           3, 2, 3, 2, 5, 3, 3, -1, 5, 2, 11, 2, 7, 5, 3, 3, 7, 2, 5, 2, 3, 7, 13, 5, 3, 2, -1);
begin
  f := 0;
  cnt := 0;
  write('Test mp_qnr ...');
  for i:=1 to NT do begin
    mp_read_decimal(a,tests[i].ns);
    j := mp_qnr(a);
    inc(cnt);
    if j<>tests[i].qnr then begin
      inc(f);
      writeln('Failure part 1, i=',i);
    end;
  end;
  for i:=1 to 100 do begin
    mp_set_int(a,i);
    j := mp_qnr(a);
    inc(cnt);
    if j<>q100[i] then begin
      inc(f);
      writeln('Failure part 2, i=',i);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_isfundamental;
const
  {Fundamental discriminants of real quadratic fields; http://oeis.org/A003658}
  fpos: array[1..60] of byte =
          (1,5,8,12,13,17,21,24,28,29,33,37,40,41,44,53,56,57,
           60,61,65,69,73,76,77,85,88,89,92,93,97,101,104,105,
           109,113,120,124,129,133,136,137,140,141,145,149,152,
           156,157,161,165,168,172,173,177,181,184,185,188,193);
const
  {Negative discriminants  Disk, http://oeis.org/A003657}
  fneg: array[1..60] of byte =
          (3,4,7,8,11,15,19,20,23,24,31,35,39,40,43,47,51,52,
           55,56,59,67,68,71,79,83,84,87,88,91,95,103,104,107,
           111,115,116,119,120,123,127,131,132,136,139,143,148,
           151,152,155,159,163,164,167,168,179,183,184,187,191);
var
  i,cnt,f: integer;
  tv: array[1..193] of boolean;
begin
  f := 0;
  cnt := 0;
  writeln('Test is_fundamental32 ...');
  fillchar(tv,sizeof(tv),ord(false));
  for i:=1 to 60 do tv[fpos[i]] := true;
  for i:=1 to 193 do begin
    inc(cnt);
    if is_fundamental32(i)<>tv[i] then begin
      writeln('failed for ',i);
      inc(f);
    end;
  end;
  fillchar(tv,sizeof(tv),ord(false));
  for i:=1 to 60 do tv[fneg[i]] := true;
  for i:=1 to 193 do begin
    inc(cnt);
    if is_fundamental32(-i)<>tv[i] then begin
      writeln('failed for ',-i);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_issquarefree;
var
  i,cnt,f: integer;
const
  tv: array[1..100] of byte = {pari: vector(100,n, issquarefree(n)}
        (1,1,1,0,1,1,1,0,0,1,1,0,1,1,1,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,
         1,1,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1,1,0,0,1,1,
         1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,1,0,1,1,1,0,1,0,1,0,1,1,1,0,1,0,0,0);
begin
  f := 0;
  cnt := 0;
  writeln('Test is_squarefree32 ...');
  for i:=1 to 100 do begin
    inc(cnt,2);
    if is_squarefree32(i) <> (tv[i]=1) then begin
      writeln('failed for i=',i);
      inc(f);
    end;
    if is_squarefree32(-i) <> (tv[i]=1) then begin
      writeln('failed for i=',-i);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_quaddisc;
var
  i,cnt,f: integer;
const
  dv: array[-50..50] of integer = ( {Pari: vector(101, n, quaddisc(n-51))}
        -8, -4, -3, -47, -184, -20, -11, -43, -168, -164, -40, -39, -152,
        -148, -4, -35, -136, -132, -8, -31, -120, -116, -7, -3, -104, -4,
        -24, -23, -88, -84, -20, -19, -8, -68, -4, -15, -56, -52, -3, -11,
        -40, -4, -8, -7, -24, -20, -4, -3, -8, -4, 0, 1, 8, 12, 1, 5, 24,
        28, 8, 1, 40, 44, 12, 13, 56, 60, 1, 17, 8, 76, 5, 21, 88, 92, 24,
        1, 104, 12, 28, 29, 120, 124, 8, 33, 136, 140, 1, 37, 152, 156,
        40, 41, 168, 172, 44, 5, 184, 188, 12, 1, 8);
begin
  f := 0;
  cnt := 0;
  writeln('Test quaddisc32 ...');
  for i:=-50 to 50 do begin
    inc(cnt);
    if quaddisc32(i) <> dv[i] then begin
      writeln('failed for i=',i);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_core32;
var
  i,cnt,f: integer;
const
  cv: array[1..77] of byte = ({http://oeis.org/A007913}
        1,2,3,1,5,6,7,2,1,10,11,3,13,14,15,1,17,2,19,5,21,
        22,23,6,1,26,3,7,29,30,31,2,33,34,35,1,37,38,39,10,
        41,42,43,11,5,46,47,3,1,2,51,13,53,6,55,14,57,58,59,
        15,61,62,7,1,65,66,67,17,69,70,71,2,73,74,3,19,77);
begin
  f := 0;
  cnt := 0;
  writeln('Test core32 ...');
  for i:=1 to 77 do begin
    inc(cnt);
    if core32(i) <> cv[i] then begin
      writeln('failed for i=',i);
      inc(f);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_rad32;
const
  rv: array[0..78] of byte = ({http://oeis.org/A007947, plus rad(0)=0}
        0,1,2,3,2,5,6,7,2,3,10,11,6,13,14,15,2,17,6,19,10,
        21,22,23,6,5,26,3,14,29,30,31,2,33,34,35,6,37,38,
        39,10,41,42,43,22,15,46,47,6,7,10,51,26,53,6,55,
        14,57,58,59,30,61,62,21,2,65,66,67,34,69,70,71,6,
        73,74,15,38,77,78);
  (* Pari: rad(n)={local(p); p=factor(n)[, 1]; prod(i=1, length(p), p[i])} *)
  ix: array[1..4] of longint = (-111111111, 123456780, -98765000, 98765400);
  iy: array[1..4] of longint = ( -37037037,  20576130,   -197530,  4938270);
var
  i,cnt,f: integer;
  r: longint;
begin
  f := 0;
  cnt := 0;
  writeln('Test rad32 ...');
  for i:=0 to 78 do begin
    inc(cnt);
    if rad32(i)<>rv[i] then begin
      inc(f);
      writeln('Fail for ',i);
    end;
  end;
  for i:=1 to 4 do begin
    inc(cnt);
    r := rad32(ix[i]);
    if r<>iy[i] then begin
      inc(f);
      writeln('Fail for ',ix[i]);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_tau32;
const
  rv: array[1..99] of byte = ({http://oeis.org/A000005}
        1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6,4,4,2,8,3,4,4,6,2,8,2,
        6,4,4,4,9,2,4,4,8,2,8,2,6,6,4,2,10,3,6,4,6,2,8,4,8,4,4,2,12,2,
        4,6,7,4,8,2,6,4,8,2,12,2,4,6,6,4,8,2,10,5,4,2,12,4,4,4,8,2,12,
        4,6,4,4,4,12,2,6,6);
  ix: array[1..4] of longint = (-9, 2147483646, 1234567890, 2095133040);
  iy: array[1..4] of integer = (3,  192, 48,  1600);
var
  i,cnt,f: integer;
  r: longint;
begin
  f := 0;
  cnt := 0;
  writeln('Test tau32 ...');
  for i:=1 to 99 do begin
    inc(cnt);
    if tau32(i)<>rv[i] then begin
      inc(f);
      writeln('Fail for ',i);
    end;
  end;
  for i:=1 to 4 do begin
    inc(cnt);
    r := tau32(ix[i]);
    if r<>iy[i] then begin
      inc(f);
      writeln('Fail for ',ix[i]);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_safeprime32;
const
  ix: array[1..5] of longint = (1, 24, 123456789, 2000000000, 2147483555);
  iy: array[1..5] of longint = (5, 47, 123457259, 2000000579, 2147483579); {Maple}
var
  i,cnt,f: integer;
  r: longint;
begin
  f := 0;
  cnt := 0;
  writeln('Test safeprime32 ...');
  for i:=1 to 5 do begin
    inc(cnt);
    r := safeprime32(ix[i]);
    if r<>iy[i] then begin
      inc(f);
      writeln('Fail for ',ix[i]);
    end;
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_valuation;
const
  {Pari: vector(50, n, valuation((47!)^29,n+1))}
  vt1: array[2..51] of word = (1218, 609, 609, 290, 609, 174, 406, 304,
         290, 116, 609, 87, 174, 290, 304, 58, 304, 58, 290, 174, 116,
         58, 406, 145, 87, 203, 174, 29, 290, 29, 243, 116, 58, 174, 304,
         29, 58, 87, 290, 29, 174, 29, 116, 290, 58, 29, 304, 87, 145, 58);
  {Pari: vector(20, n, valuation((400!)^3,n+256))}
  vt2: array[1..20] of word = (3, 27, 30, 96, 39, 9, 3, 117, 21, 66,
         12, 15, 3, 196, 3, 72, 96, 6, 117, 51);
var
  i,j,cnt,f: integer;
  v: longint;
begin
  writeln('Test mp_val/rem ');
  cnt := 0;
  f := 0;
  mp_fact(47,a);
  mp_expt_w(a,29,a);
  for i:=2 to 51 do begin
    write('.');
    inc(cnt);
    v := mp_val(a,i);
    if v<>vt1[i] then begin
      inc(f);
      writeln('mp_val failed for i=',i);
    end;
    inc(cnt);
    mp_valrem(a,-i,b,v);
    if v<>vt1[i] then begin
      inc(f);
      writeln('mp_valrem: wrong v for i=',i);
    end
    else begin
      {test a = b*i^v}
      mp_set_pow(c,-i,v);
      mp_mul(c,b,d);
      if mp_is_ne(a,d) then begin
        inc(f);
        writeln('mp_valrem failed for i=',i);
      end;
    end;
  end;
  {check non-digit version for MP_16BIT}
  mp_fact(400,a);
  mp_expt_w(a,3,a);
  for i:=1 to 20 do begin
    write('.');
    inc(cnt);
    j := i+256;
    v := mp_val(a,j);
    if v<>vt2[i] then begin
      inc(f);
      writeln('mp_val failed for i=',i);
    end;
    inc(cnt);
    mp_valrem(a,-j,b,v);
    if v<>vt2[i] then begin
      inc(f);
      writeln('mp_valrem: wrong v for i=',i);
    end
    else begin
      {test a = b*i^v}
      mp_set_pow(c,-j,v);
      mp_mul(c,b,d);
      if mp_is_ne(a,d) then begin
        inc(f);
        writeln('mp_valrem failed for i=',i);
      end;
    end;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_nroot_modp;
var
  err,f: integer;

  procedure setabc(ia,ib,ic: longint);
  begin
    mp_set_int(a,ia);
    mp_set_int(b,ib);
    mp_set_int(c,ic);
  end;

  procedure check(n,ie: longint);
  begin
    if err=0 then begin
      if mp_cmp_int(e,ie)<>0 then begin
        inc(f);
        writeln(' Test ',n,'  Diff: ', ie, ' <> ', mp_decimal(e));
      end;
    end
    else begin
      writeln(' Test ',n,'  error=',err);
      inc(f);
    end;
  end;

  procedure checkerr(n,ie: longint);
  begin
    if err<>ie then begin
      inc(f);
      writeln(' Test ',n,'  Unexpected error:e ', err, ' <> ', ie);
    end;
  end;

begin
  writeln('Test mp_is_power_modpk and s_mp_nroot_modp/pk/pq');
  f := 0;
  {632^101 = 123 mod 997}
  setabc(123,101,997);
  s_mp_nroot_modp(a,b,c,e,err);
  check(1,632);

  {295^(-101) = 123 mod 997}
  setabc(123,-101,997);
  s_mp_nroot_modp(a,b,c,e,err);
  check(2,295);

  {x^10 = 9 mod 41 solvable, e.g. x=3}
  setabc(9,10,41);
  s_mp_nroot_modp(a,b,c,e,err);
  checkerr(3,3);
  if not mp_is_power_modpk(a,b,c,1) then begin
    inc(f);
    writeln('Test 4 failed');
  end;

  {23563^101 = 123 mod 257*997}
  setabc(123,101,257);
  mp_set(d,997);
  s_mp_nroot_modpq(a,b,c,d,e,err);
  check(5,23563);

  {48151^(-101) = 123 mod 257*997}
  mp_set_int(b,-101);
  s_mp_nroot_modpq(a,b,c,d,e,err);
  check(6,48151);

  {14^7=4 mod 25}
  setabc(4,7,5);
  s_mp_nroot_modpk(a,b,c,2,e,err);
  check(7,14);

  {12^(-7) = 22 mod 25}
  setabc(22,-7,5);
  s_mp_nroot_modpk(a,b,c,2,e,err);
  check(8,12);

  {x^5 = 7 mod 5^2 solvable, e.g. x=2}
  {gcd(n,phi(p^k)) = gcd(5,20) = 5}
  setabc(7,5,5);
  s_mp_nroot_modpk(a,b,c,2,e,err);
  checkerr(9,2);
  if not mp_is_power_modpk(a,b,c,2) then begin
    inc(f);
    writeln('Test 10 failed');
  end;

  {x^890123 = 134567 mod 1234567891^7 with}
  {x = 1736114324929414064898236186182891179804698965437530215311812105}
  setabc(134567, 890123, 1234567891);
  mp_read_decimal(d,'1736114324929414064898236186182891179804698965437530215311812105');
  s_mp_nroot_modpk(a,b,c,7,e,err);
  if err=0 then begin
    if mp_cmp(e,d)<>0 then begin
      inc(f);
      writeln(' Test 11 - Diff: ', mp_decimal(d), ' <> ', mp_decimal(e));
    end;
  end
  else begin
    inc(f);
    writeln(' Test 11 - error=',err);
  end;

  {x^5 = 27 mod 41 solvable, e.g. x=14}
  {gcd(b,c-1)) = gcd(5,40) = 5}
  setabc(27,5,41);
  s_mp_nroot_modp(a,b,c,e,err);
  check(12,14);

  {x^11 = 58 mod 199 with x=43 and others}
  setabc(58,11,199);
  s_mp_nroot_modp(a,b,c,e,err);
  check(13,43);

  setabc(57,11,199);
  if mp_is_power_modpk(a,b,c,1) then begin
    inc(f);
    writeln('Test 14 failed');
  end;

  writeln('No. of checks: ',14,', failed: ',f);
  inc(totalfailed,f);
  writeln;

end;


{---------------------------------------------------------------------------}
procedure test_primepower32;
var
  n,i,k,f: longint;
  r: integer;
  bl,b3: boolean;
{$ifdef BIT16}
const
  CNT = 100000;
{$else}
const
  CNT = 150000;
{$endif}
begin
  writeln('Test is_primepower32 (compare with mp_is_primepower)');
  f := 0;
  for n :=1 to CNT do begin
    if n and $7FF = 0 then write('.');
    mp_set_int(a,n);
    bl := mp_is_primepower(a,b,i);
    r  := is_primepower32(n,k);
    b3 := r>1;
    if bl<>b3 then begin
      if f=0 then writeln(' fail for: ',n);
      inc(f);
    end;
  end;
  writeln;
  writeln('No. of checks: ',CNT,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_safeprime_ex;
var
  p,g: longint;
  i,f: integer;
{$ifdef BIT16}
const
  CNT = 500;
  DELTA = 10000;
{$else}
const
  CNT = 2000;
  DELTA = 10000;
{$endif}
begin
  writeln('Test mp_safeprime_ex');
  p := 1;
  f := 0;
  mp_zero(a);
  for i:=1 to CNT do begin
    if i and $1F = 0 then write('.');
    p := safeprime32(p+DELTA);
    g := primroot32(p);
    mp_add_d(a,DELTA,a);
    mp_safeprime_ex(a,b,true);
    mp_set_int(c,p);
    mp_set_int(d,g);
    if mp_is_ne(a,c) or mp_is_ne(b,d) then begin
      writeln(p:8, g:8);
    end;
  end;
  writeln;
  writeln('No. of checks: ',CNT,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_dlog32;
var
  i,f,cnt: integer;
  n,g,b,x,y: longint;
begin
  {Note: More tests in t_dlog32.pas in supptest.zip}
  writeln('Test dlog32');
  f := 0;
  cnt := 0;

  n := 107;
  g := 2;
  {2 is a primitive root mod 107, there should be no failure here}
  for b:=1 to n-1 do begin
    x := dlog32(g,b,n);
    y := exptmod32(g,x,n);
    inc(cnt);
    if (x=-1) or (y<>b) then begin
      inc(f);
      writeln(' Test ', cnt,' failed', b:4, x:12, y:12, y=b:10);
    end;
  end;

  {3 is no primitive root mod 107, but 3^15 = 100 mod 107}
  {and 3^x = 98 mod 107 has no solution}
  inc(cnt);
  if dlog32(3,100,n)<>15 then begin
    inc(f);
    writeln('Test 110 failed');
  end;
  inc(cnt);
  if dlog32(3, 98,n)<>-1 then begin
    inc(f);
    writeln('Test 111 failed');
  end
  else write('.');

  {HAC example: 2 is no primitive root mod 383, but 2^110 = 228 mod 383}
  n := 383;
  g := 2;
  b := 228;
  x := dlog32(g,b,n);
  y := exptmod32(g,x,n);
  inc(cnt);
  if y<>b then begin
    inc(f);
    writeln('Test 112 failed:', x:12, y:12, y=b:10);
  end
  else write('.');

  n := 1048703;
  g := 5;
  b := 228;
  x := dlog32(g,b,n);
  y := exptmod32(g,x,n);
  inc(cnt);
  if y<>b then begin
    inc(f);
    writeln('Test 113 failed:', x:12, y:12, y=b:10);
  end
  else write('.');

  {must fail g=-1!}
  n := 1048703;
  g := n-1;
  b := 228;
  x := dlog32(g,b,n);
  inc(cnt);
  if x<>-1 then begin
    inc(f);
    writeln('Test 114 failed: Spurious solution found');
  end
  else write('.');

  n := 536871263;
  g := 5;
  b := 228;
  x := dlog32(g,b,n);
  y := exptmod32(g,x,n);
  inc(cnt);
  if y<>b then begin
    inc(f);
    writeln('Test 115 failed:', x:12, y:12, y=b:10);
  end
  else write('.');

  n := 1073742623;
  g := 5;
  b := 228;
  x := dlog32(g,b,n);
  y := exptmod32(g,x,n);
  inc(cnt);
  if y<>b then begin
    inc(f);
    writeln('Test 116 failed:', x:12, y:12, y=b:10);
  end
  else write('.');

  n := 1073732879;
  g := 7;
  b := 228;
  x := dlog32(g,b,n);
  y := exptmod32(g,x,n);
  inc(cnt);
  if y<>b then begin
    inc(f);
    writeln('Test 117 failed:', x:12, y:12, y=b:10);
  end
  else write('.');

  {Largest safe prime < 2^31}
  n := 2147483579;
  g := 2;
  for b:=2 to 36 do begin
    x := dlog32(g,b,n);
    y := exptmod32(g,x,n);
    inc(cnt);
    if y<>b then begin
      inc(f);
      writeln('Test ',120+b,' failed:', x:12, y:12, y=b:10);
    end
    else write('.');
  end;

  {Some random powers of smallest primitive root for random primes}
  n := 12345678;
  for i:=0 to 20 do begin
    n := nextprime32(n + random($FFFF));
    g := primroot32(n);
    b := exptmod32(g, 100 + random(10000),n);
    x := dlog32(g,b,n);
    y := exptmod32(g,x,n);
    inc(cnt);
    if y<>b then begin
      inc(f);
      writeln('Test ',150+i,' failed:', x:12, y:12, y=b:10);
    end
    else write('.');
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_is_primeroot32;
var
  f,cnt: integer;
  k,n,g,pc: longint;
{$ifdef BIT16}
const
  NMAX = 400;
{$else}
const
  NMAX = 1000;
{$endif}
begin
  writeln('Test is_primeroot32');
  f := 0;
  cnt := 0;
  for n:=1 to NMAX do begin
    {n has primitive roots if n = 2, 4, p^k, 2p^k}
    if n and 7 = 0 then continue;
    if (is_primepower32(n,g)>0) or ((n and 1 = 0) and (is_primepower32(n div 2,g)>0)) then begin
      inc(cnt);
      if cnt and 3 = 0 then write('.');
      pc := 0;
      for k:=1 to n do begin
        if is_primroot32(k,n) then inc(pc);
      end;
      k := EulerPhi32(EulerPhi32(n));
      if pc<>k then begin
        inc(f);
        writeln(' Failure for n=',n);
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

label
  done;
begin
  writeln('Test of MPArith V', MP_VERSION, ' [main test]  (c) 2004-2018 W.Ehrhardt ');
  writeln('MAXDigits = ',MAXDigits, ',  MP_MAXBIT = ',MP_MAXBIT, ',  MAXRadix = ', MAXRadix);
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  writeln;

  mp_init5(a,b,c,d,e);
  mp_init3(s,t,z);
  mp_rsa_init_private(prk);

  totalfailed := 0;

{
  reverse_check;
  goto done;
}
  testradixconv;
  if totalfailed<>0 then begin
    writeln('Failures in testradixconv! Program aborted!');
    goto done;
  end;

  {$ifndef BIT16}
    astr_test;
  {$endif}

  mp_show_plus := false;
  mp_uppercase := false;

  test_checksum;

  TestPrimes16Index;
  test_pi_lphi;
  test_prime32;
  gcd32_test;
  xgcd32_test;
  test_carmichael;
  test_EulerPhi;
  test_order32;
  test_primroot;
  test_moebius;
  test_sigmak;
  test_isfundamental;
  test_issquarefree;
  test_core32;
  test_quaddisc;
  test_rad32;
  test_tau32;
  test_safeprime32;
  test_primepower32;
  test_is_primeroot32;
  test_dlog32;

  IsLong_Test;
  add_test;
  sub_test;
  mul_test1;
  mul_test2;
  sqr_test;
  mod_test;
  div_test;
  expt_test;
  reduce_test;
  gr_mod_test;
  invmod_test;
  emod_test;
  emod2_test;
  red2k_test;
  egcd_test;
  kron_test;
  jac_test;
  lucas_test;
  coshmult_test;
  lucfib_test;
  test_lucasuv;
  test_catalan;
  pell1_test;
  pell4_test;
  test_rqffu;
  powerd_test;
  cornacchia_test;
  test_4sq;
  bool_test;
  nroottest;
  facttest;
  swing_test;
  binomial_test;
  test_poch_perm;
  probprime_test;
  test_rnr;
  rsa_test1;
  rsa_test2;
  rsa_test3;
  sqrtmod_test;
  sqrtmod2k_test;
  sqrtmodpq_test;
  cbrtmod_test;
  cbrtmodpq_test;
  test_nroot_modp;;
  test_squad_mod;
  test_qnr;
  test_valuation;
  test_safeprime_ex;
  test_is_carmichael;
  digitsum_test;
  reverse_check;

  if paramcount=0 then pfdu_test;

done:

  mp_rsa_clear_private(prk);
  mp_clear5(a,b,c,d,e);
  mp_clear3(s,t,z);

  if totalfailed<>0 then begin
    writeln('************************');
    writeln('****** Total failed: '#7#7#7, totalfailed);
    writeln('************************');
  end
  else writeln('All test OK: no failures');

  writeln;

  {$ifdef MPC_Diagnostic}
    mp_dump_meminfo;
    mp_dump_diagctr;
  {$endif}
end.

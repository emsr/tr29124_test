unit mp_prime;

{Basic 16/32 bit prime numbers and functions}

interface

{$i STD.INC}

{$ifdef VER70}
{$A+,N+}
{$endif}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  Basic 16/32 bit prime numbers and functions

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  mp_prm16.inc, mp_pbits.inc

 MEMORY USAGE  :  4 KB pbits, 13KB Primes16
                  heap for prime sieve routines

 DISPLAY MODE  :  ---

 REFERENCES    :  [3] Knuth, D.E.: The Art of computer programming. Vol 2
                      Seminumerical Algorithms, 3rd ed., 1998
                  [4] Forster, O.: Algorithmische Zahlentheorie, 1996
                  [5] (HAC) Menezes,A., von Oorschot,P., Vanstone, S: Handbook of
                      Applied Cryptography, 1996, www.cacr.math.uwaterloo.ca/hac
                  [7] P. Ribenboim: The New Book of Prime Number Records, 3rd ed., 1995.
                  [8] Marcel Martin: NX - Numerics library of multiprecision
                      numbers for Delphi and Free Pascal, 2006-2009
                      www.ellipsa.eu/public/nx/index.html
                 [10] R.Crandall, C.Pomerance: Prime Numbers, A Computational
                      Perspective, 2nd ed., 2005
                 [12] PARI/GP at http://pari.math.u-bordeaux.fr/
                 [24] H. Cohen, A Course in Computational Algebraic Number Theory
                      4th printing, 2000
                 [39] D.H. Lehmer: On the exact number of primes less than a given limit,
                      Illinois Journal of Mathematics, vol. 3, pp. 381-388, 1959;
                      available from http://projecteuclid.org/euclid.ijm/1255455259
                 [40] H. Riesel, Prime Numbers and Computer Methods for Factorization,
                      Vol. 126 of Progress in Mathematics, Boston, 2nd ed. 1994.
                      Paperback reprint 2012 in Modern Birkh„user Classics Series.



 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.9.00   26.12.08  W.Ehrhardt  Initial version: IsPrime16/32,is_spsp32,is_spsp32A from mp_base
 1.9.01   29.12.08  we          prime_sieve routines
 1.9.02   30.12.08  we          prime_sieve with mp_getmem, prime_sieve_reset
 1.9.03   04.01.09  we          Primes16 array as global data, Primes16Index
 1.9.04   04.01.09  we          32 bit first/next prime routines from mp_numth

 1.12.00  05.07.09  we          Removed prime residue classes mod 30

 1.20.00  14.01.12  we          BIT64: _spsp32, is_spsp32A, IsPrime32, modnpd2
 1.20.01  21.01.12  we          IsPrime16 inline if supported

 1.21.00  08.07.12  we          Assert psieve<>nil and paux<>nil in prime_sieve_next
 1.21.01  14.07.12  we          lsumphi32 and primepi32
 1.21.02  15.07.12  we          improved Primes16Index
 1.21.03  23.07.12  we          TFactors32 and PrimeFactor32
 1.21.04  23.07.12  we          EulerPhi32 and Carmichael32
 1.21.05  24.07.12  we          Moebius32
 1.21.06  25.07.12  we          primroot32
 1.21.07  27.07.12  we          order32

 1.22.00  06.08.12  we          SIEVE_MAXPRIME = MaxLongint if not MPC_SmallSieve
 1.22.01  07.08.12  we          prime32
 1.22.02  03.09.12  we          is_fundamental32
 1.22.03  04.09.12  we          is_squarefree32, improved Moebius32
 1.22.04  05.09.12  we          primepi32 uses icbrt32
 1.22.05  05.09.12  we          core32, quaddisc32, rad32

 1.24.00  04.01.13  we          BIT64 changed to PurePascal

 1.29.00  14.07.14  we          tau32
 1.29.01  15.07.14  we          safeprime32

 1.30.00  28.09.14  we          is_primepower32
 1.30.01  30.09.14  we          dlog32
 1.30.02  01.10.14  we          dlog32_ex, expanded and tuned version
 1.30.03  03.10.14  we          is_primroot32
 1.30.04  04.10.14  we          dlog32/ex without safe prime requirement

 1.33.00  08.09.16  we          is_Carmichael32

 1.38.00  25.06.18  we          prime32 uses doubles

**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2018 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)


{#Z+}
{---------------------------------------------------------------------------}
{------------------------ Basic prime functions ----------------------------}
{---------------------------------------------------------------------------}
{#Z-}

const
  NumPrimes16 = 6542;                      {Number of 16 bit primes}

const
  Primes16: array[1..NumPrimes16] of word  {array of all 16 bit primes}
               = {#Z+}({$i mp_prm16.inc}){#Z-};

const
  pbits16: array[0..4095] of byte = {#Z+}({$i mp_pbits.inc}){#Z-};
           {Bit array for 16 bit primes. Only odd integers are represented}
  pmask16: array[0..15] of byte = (0,1,0,2,0,4,0,8,0,16,0,32,0,64,0,128);
           {Mask array: odd i is prime if pbits16[i shr 4] and pmask16[i and $f] <> 0)}


function IsPrime16(N: word): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Test if N is prime}

function Primes16Index(n: word): integer;
  {-Return index of largest prime <= n in Primes16 array; 1 if n<=2}
  { Since Primes16[1]=2, this is identical to primepi(n) for n>=2. }

function IsPrime32(N: longint): boolean;
  {-Test if longint (DWORD) N is prime}

function is_primepower32(n: longint; var p: longint): integer;
  {-Test if n is a prime power: return 0 if not, n = p^result otherwise.}
  { Note: contrary to mp_is_primepower a prime n is returned as n = n^1!}

function lsumphi32(x: longint; a: integer): longint;
  {-Return the partial sieve function Phi(x,a), the number of integers in}
  { (0,x] which are co-prime to the first a primes, aka 'Legendre's sum'}

function primepi32(x: longint):  longint;
  {-Return the prime counting function pi(x) using Lehmer's formula}

function is_spsp32(N, A: longint): boolean;
  {-Strong probable prime test for N with base A, calls is_spsp32A}

function is_spsp32A(N: longint; const bases: array of longint): boolean;
  {-Strong probable prime test for N with a number of bases. }
  { Return true if N is a probable prime for all bases, or if}
  { N=2. Negative numbers are treated as unsigned (=cardinal)}

{#Z+}
{---------------------------------------------------------------------------}
{--------------------- 32-bit number-theoretic functions -------------------}
{---------------------------------------------------------------------------}
{#Z-}

type
  {Record for prime factorization, 10 primes are enough for 32 bit numbers,}
  {because numbers with 10 different prime factors are >= 6469693230 > 2^32}
  {number = product(primes[i]^pexpo[i], i=1..pcount), if n > 1}
  TFactors32 = record
                 number: longint;                 {original number}
                 primes: array[1..10] of longint; {prime factors in ascending order}
                 pexpo : array[1..10] of byte;    {exponents of prime factors}
                 pcount: integer;                 {number of different primes}
               end;


function  Carmichael32(n: longint): longint;
  {-Return the Carmichael function lambda(|n|), lambda(0)=0. For n > 0 this}
  { is the least k such that x^k = 1 (mod n) for all x with gcd(x,n) = 1.  }

function  core32(n: longint): longint;
  {-Return the squarefree part of n, n<>0, i.e. the unique squarefree integer c with n=c*f^2}

function  dlog32(a,b,p: longint): longint;
  {-Compute the discrete log_a(b) mod p using Pollard's rho algorithm: i.e.}
  { solve a^x = b mod p, with p prime, a > 1, b > 0; return -1 for failure.}
  { If a is no primitive root mod p, solutions may not exist or be unique. }

function  dlog32_ex(a,b,p,JMAX: longint): longint;
  {-Compute the discrete log_a(b) mod p using Pollard's rho  }
  { algorithm; p prime, a > 1, b > 0; return -1 for failure. }
  { Expanded version with variable trial parameter JMAX >= 0.}

function  EulerPhi32(n: longint): longint;
  {-Return Euler's totient function phi(|n|), phi(0)=0. For n > 0 }
  { this the number of positive integers k <= n with gcd(n,k) = 1.}

function  is_Carmichael32(n: longint): boolean;
  {-Test if |n| is a Carmichael number}

function  is_fundamental32(d: longint): boolean;
  {-Return true, if d is a fundamental discriminant (either d=1 mod 4 and }
  { squarefree, or d=0 mod 4, d/4 = 2,3 mod and squarefree), false if not.}

function  is_primroot32(const g,n: longint): boolean;
  {-Test if g is primitive root mod n}

function  is_squarefree32(n: longint): boolean;
  {-Return true if n is squarefree, i.e. not divisible by a square > 1}

function  Moebius32(n: longint): integer;
  {-Return the Moebius function mu(abs(n)), mu(0)=0, mu(1)=1. mu(n)=(-1)^k }
  { if n > 1 is the product of k different primes; otherwise mu(n)=0.      }

function  order32(n,m: longint): longint;
  {-Return the order of n mod m, m > 1, i.e. the smallest integer}
  { e with n^e = 1 mod m; if gcd(n,m) <> 1 the result is 0.}

function  prime32(k: longint): longint;
  {-Return the kth prime if 1 <= k <= 105097565, 0 otherwise}

procedure PrimeFactor32(n: longint; var FR: TFactors32);
  {-Return the prime factorization of n > 1, FR.pcount=0 if n < 2}

function  primroot32(n: longint): longint;
  {-Compute the smallest primitive root mod n, 0 if n does not have a prim.root}

function  quaddisc32(n: longint): longint;
  {-Return the discriminant of the quadratic field Q(sqrt(n))}

function  rad32(n: longint): longint;
  {-Return the radical rad(n) of n = product of the distinct prime factors of n.}

function  tau32(n: longint): integer;
  {-Return the number of positive divisors of n (including 1 and n)}


{#Z+}
{---------------------------------------------------------------------------}
{------------------------ Prime sieve functions ----------------------------}
{---------------------------------------------------------------------------}
{#Z-}

type
  TPrimeContext = record             {Context for FindFirst/NextPrime32}
                    prime: longint;  {last found prime}
                    next : longint;  {next value to test}
                    index: integer;  {index in residue class table}
                  end;

procedure FindFirstPrime32(n: longint; var ctx: TPrimeContext);
  {-Find first prime >= n and initialize ctx for FindNextPrime32}

procedure FindNextPrime32(var ctx: TPrimeContext);
  {-Find next 32 bit prime (DWORD interpretation, see note), success if ctx.prime<>0}

function  nextprime32(n: longint): longint;
  {-Next 32 bit prime >= n (DWORD interpretation, see note)}

procedure nextprime32_array(n,cmax: longint; var a: array of longint);
  {-Fill an array with the next 32 bit primes >= n (DWORD interpretation).}
  { If a[i0]=4294967291 then a[i]=0 for i>i0. If cmax<=0 fill complete array.}

function  prevprime32(n: longint): longint;
  {-Previous 32 bit prime <= n, prevprime32(0)=0, (DWORD interpretation)}

function  safeprime32(n: longint): longint;
  {-Return the next safe prime p >= n, i.e. p and (p-1)/2 prime; 0 if n > 2147483579}

{$ifdef MPC_SmallSieve}
const
  SIEVE_MAXPRIME = 1908867043; {maximum prime that will be generated, see note}
  SIEVE_PRIMES   = 4550;       {= primepi(sqrt(SIEVE_MAXPRIME))-1}
  SIEVE_BLOCK    = 32768;      {sieve size, may be set to 16384 for BIT16 if memory is low}

type
  TAux   = array[0..SIEVE_PRIMES-1] of word;   {prime offsets in sieve}
  TFlags = array[0..SIEVE_BLOCK-1] of boolean; {prime flags in sieve}
  TSieve = record                              {prime sieve context}
             paux    : ^TAux;                  {pointer to offsets }
             psieve  : ^TFlags;                {pointer to flags   }
             curr_blk: longint;                {current sieve block}
             curr_off: longint;                {offset in curr_blk }
           end;
{$else}
const
  SIEVE_MAXPRIME = MaxLongint; {maximum prime that will be generated, see note}
  SIEVE_PRIMES   = 4792;       {= primepi(sqrt(SIEVE_MAXPRIME))-1}
  SIEVE_BLOCK    = 32768;      {sieve size, may be set to 16384 for BIT16 if memory is low}

type
  TAux   = array[0..SIEVE_PRIMES-1] of longint;{prime offsets in sieve}
  TFlags = array[0..SIEVE_BLOCK-1] of boolean; {prime flags in sieve}
  TSieve = record                              {prime sieve context}
             paux    : ^TAux;                  {pointer to offsets }
             psieve  : ^TFlags;                {pointer to flags   }
             curr_blk: longint;                {current sieve block}
             curr_off: longint;                {offset in curr_blk }
           end;
{$endif}


procedure prime_sieve_clear(var sieve: TSieve);
  {-Release memory allocated by prime_sieve_init}

procedure prime_sieve_init(var sieve: TSieve; first_prime: longint);
  {-Allocate/initialize sieve to return primes >= first_prime; first_prime <= SIEVE_MAXPRIME}

function  prime_sieve_next(var sieve: TSieve): longint;
  {-Return next prime from sieve, 1 if done}

procedure prime_sieve_reset(var sieve: TSieve; first_prime: longint);
  {-Initialize already allocated sieve to return primes >= first_prime; first_prime <= SIEVE_MAXPRIME}


{#Z+}
{Tables for prime residue classes mod 210, calculated with t_rcnp.pas}
{Must be interfaced because they are used by mp_nextprime/mp_prevprime.}
{NPRC_OddIdx has the number of the prime residue classes of an odd value }
{n mod NPRC_MOD, or -1 if there is no corresponding prime residue class. }
{NPRC_Diff contains the differences: NPRC_Diff[i] = prctab[i+1]-prctab[i]}

{prime residue classes mod 210. Cnt=48=1*2*4*6, numbered from 0 to 47}
  { 1,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,
   53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103,
  107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157,
  163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209}

const
  NPRC_MOD = 210;
  NPRC_NRC = 48;
const
  NPRC_OddIdx:  array[0..(NPRC_MOD div 2)-1] of shortint =
                ( 0, -1, -1, -1, -1,  1,  2, -1,  3,  4, -1,  5, -1, -1,  6,
                  7, -1, -1,  8, -1,  9, 10, -1, 11, -1, -1, 12, -1, -1, 13,
                 14, -1, -1, 15, -1, 16, 17, -1, -1, 18, -1, 19, -1, -1, 20,
                 -1, -1, -1, 21, -1, 22, 23, -1, 24, 25, -1, 26, -1, -1, -1,
                 27, -1, -1, 28, -1, 29, -1, -1, 30, 31, -1, 32, -1, -1, 33,
                 34, -1, -1, 35, -1, -1, 36, -1, 37, 38, -1, 39, -1, -1, 40,
                 41, -1, -1, 42, -1, 43, 44, -1, 45, 46, -1, -1, -1, -1, 47);
const
  NPRC_Diff:    array[0..NPRC_NRC-1] of mp_digit =
                (10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4,
                  2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4,
                  6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2,10, 2);
{#Z-}


implementation


uses
  {$ifdef MPC_PurePascal} BTypes,{$endif}
  mp_base, mp_prng;


{---------------------------------------------------------------------------}
{------------------------ Basic prime functions ----------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function IsPrime16(N: word): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Test if N is prime}
begin
  IsPrime16 := (N=2) or (pbits16[N shr 4] and pmask16[N and $0F] <> 0);
end;


{---------------------------------------------------------------------------}
function Primes16Index(n: word): integer;
  {-Return index of largest prime <= n in Primes16 array; 1 if n<=2}
  { Since Primes16[1]=2, this is identical to primepi(n) for n>=2.  }
var
  i,j,k,m: integer;
  p: word;
const
  js: array[0..32] of integer = (1,309,564,801,1028,1254,1469,1681,1900,2110,
                                 2312,2517,2725,2918,3124,3314,3512,3716,3908,
                                 4098,4288,4495,4678,4858,5051,5239,5432,5616,
                                 5814,6003,6179,6363,6541);
begin
  {Binary search in Primes16 array}
  Primes16Index := 1;  {Avoid FPC warning}
  if n<=2 then exit
  else if n>=65521 then Primes16Index := NumPrimes16
  else begin
    k := n shr 11;
    {Get initial lower and upper bounds better than 1 and NumPrimes16}
    i := js[k];
    j := js[k+1];
    m := 0;
    repeat
      k := (i+j) shr 1;
      if k=m then break;
      p := Primes16[k];
      if n < p then j := k-1
      else if n > p then i := k+1
      else break;
      m := k;
    until false;
    Primes16Index := k;
  end;
end;


{$ifdef BIT16}
  {$undef  USEA5}  {smaller but slower lsumphi32}
{$else}
  {$define USEA5}  {faster  but larger lsumphi32}
{$endif}


{---------------------------------------------------------------------------}
function lsumphi32(x: longint; a: integer): longint;
  {-Return the partial sieve function Phi(x,a), the number of integers in}
  { (0,x] which are co-prime to the first a primes, aka 'Legendre's sum'}
var
  sum: longint; {accumulates contributions from recursive calls}

{$ifdef USEA5}
const
  BaseA = 5;
  ProdA = 2310;
  PrA_2 = 1155;
  PrA_4 = 578;
  PhiPA = 480;
const
  {compressed table of Phi(x mod 2310, 5), TabPhiA[i] = lphi(2*i,5)}
  TabPhiA: array[0..PrA_4] of byte = (
             0, 1,1,1,1,1,1,2,2,3,4,4,5,5,5,6,7,7,7,8,8,9,10,10,11,11,11,
             12,12,12,13,14,14,14,15,15,16,17,17,17,18,18,19,19,19,20,20,
             20,20,21,21,22,23,23,24,25,25,26,26,26,26,26,26,26,27,27,28,
             28,28,29,30,30,30,30,30,31,32,32,32,33,33,33,34,34,35,36,36,
             37,37,37,38,39,39,39,39,39,40,41,41,42,43,43,43,43,43,43,44,
             44,44,44,44,45,46,46,47,48,48,49,49,49,50,51,51,51,52,52,53,
             53,53,54,54,54,55,55,55,56,57,57,57,58,58,59,60,60,60,61,61,
             62,62,62,63,63,63,63,64,64,65,66,66,67,67,67,68,68,68,68,69,
             69,69,70,70,70,70,70,71,72,72,73,73,73,74,75,75,75,76,76,76,
             77,77,78,79,79,80,80,80,81,82,82,82,83,83,84,85,85,85,86,86,
             86,86,86,87,88,88,88,88,88,89,90,90,91,92,92,93,93,93,94,94,
             94,94,95,95,96,97,97,98,98,98,98,98,98, 99, 100,100,100,101,
             101,102,103,103,103,104,104,105,105,105,106,106,106,106,106,
             106,107,108,108,109,110,110,111,111,111,111,112,112,112,113,
             113,114,114,114,115,116,116,117,117,117,118,119,119,119,120,
             120,120,120,120,121,122,122,123,123,123,124,125,125,125,126,
             126,127,128,128,129,130,130,130,130,130,131,132,132,132,132,
             132,133,134,134,135,135,135,136,136,136,137,138,138,138,139,
             139,139,140,140,141,141,141,142,142,142,143,144,144,144,145,
             145,146,147,147,147,148,148,149,149,149,150,150,150,150,151,
             151,152,153,153,153,154,154,155,155,155,155,156,156,156,157,
             157,158,158,158,159,160,160,161,161,161,162,162,162,162,163,
             163,163,164,164,165,166,166,166,166,166,167,168,168,168,169,
             169,170,171,171,172,173,173,173,173,173,174,175,175,175,175,
             175,176,177,177,178,179,179,180,180,180,180,181,181,181,182,
             182,183,184,184,185,185,185,186,186,186,187,188,188,188,189,
             189,190,190,190,190,191,191,192,192,192,193,193,193,193,194,
             194,195,196,196,197,198,198,199,199,199,199,200,200,200,201,
             201,202,202,202,203,203,203,204,204,204,205,206,206,206,207,
             207,207,208,208,209,210,210,211,211,211,212,213,213,213,214,
             214,215,216,216,217,218,218,218,218,218,219,220,220,220,220,
             220,221,222,222,222,223,223,224,224,224,225,226,226,226,227,
             227,228,229,229,230,230,230,231,231,231,232,232,232,232,233,
             233,234,235,235,235,236,236,236,236,236,237,237,237,237,238,
             238,239,240,240);
{$else}
const
  BaseA = 4;
  ProdA = 210;
  PrA_2 = 105;
  PrA_4 = 53;
  PhiPA = 48;
const
  {compressed table of Phi(x mod 210, 4), TabPhiA[i] = lphi(2*i,4)}
  TabPhiA: array[0..PrA_4] of shortint = (
             0,1,1,1,1,1,2,3,3,4,5,5,6,6,6,7,8,8,8,9,9,10,11,11,
             12,12,12,13,13,13,14,15,15,15,16,16,17,18,18,18,19,
             19,20,20,20,21,21,21,21,22,22,23,24,24);
{$endif}

  procedure phirec(x: longint; a, sign: integer);
    {-Internal function for Phi(x,a), uses the accumulator variable sum}
  begin
    {This is an optimized Pascal implementation of the C snippet given in}
    {Table II of T. Oliveira e Silva, Computing pi(x): the combinatorial }
    {method, Revista do DETUA, vol. 4, no. 6, pp. 759-768, March 2006.   }
    {Available online from http://www.ieeta.pt/~tos/bib/5.4.pdf          }
    repeat
      {T. Oliveira e Silva uses recursion downto a=0 in the demo, I stop }
      {at a=5 or a=4, and use a second function for smaller values of a. }
      if a=BaseA then begin
        {Compute Phi(x,BaseA) using the division and symmetry properties }
        {of Phi(x,a), see Lehmer[39], p.384 or Riesel[40], p.14-17       }
        a := x mod ProdA;
        x := (x div ProdA)*PhiPA;
        if a<PrA_2 then inc(x, integer(TabPhiA[succ(a) shr 1]))
        else inc(x, PhiPA - integer(TabPhiA[(ProdA-a) shr 1]));
        if sign<0 then dec(sum,x) else inc(sum,x);
        exit;
      end
      else if x < primes16[a+1] then begin
        inc(sum,sign);
        exit;
      end
      else begin
        {Use recursive property of Phi}
        dec(a);
        phirec(x div primes16[a+1], a, -sign);
      end;
    until false;
  end;

  function slowphi(x: longint; a: integer): longint;
    {-Internal function for weird cases when Phi is called with a<BaseA}
  begin
    {This is just the recursive property of Phi}
    if a<=0 then slowphi := x
    else slowphi := slowphi(x, a-1) - slowphi(x div Primes16[a], a-1);
  end;

begin
  if x<0 then lsumphi32 := -lsumphi32(-x,a)
  else if a<BaseA then lsumphi32 := slowphi(x,a)
  else begin
    sum := 0;
    phirec(x,a,1);
    lsumphi32 := sum;
  end;
end;


{---------------------------------------------------------------------------}
function primepi32(x: longint):  longint;
  {-Return the prime counting function pi(x) using Lehmer's formula}
var
  sum,z: longint;
  a,b,c,i,j,k: integer;
begin
  if x <= $FFFF then begin
    if x < 2 then primepi32 := 0
    else primepi32 := Primes16Index(word(x));
    exit;
  end;

  z := isqrt32(x);
  b := Primes16Index(z);                         {b = pi(x^(1/2) <= 4792}
  c := Primes16Index(icbrt32(x));                {c = pi(x^(1/3) <=  209}
  a := Primes16Index(isqrt32(z));                {a = pi(x^(1/4) <=   47}

  sum := longint(b+a-2)*(b-a+1) shr 1;
  for i:=a+1 to c do begin
    z := x div Primes16[i];
    sum := sum - primepi32(z);
    k := Primes16Index(isqrt32(z));
    for j:=i to k do begin
      {z div Primes16[j] <= x / Primes16[a+1]^2 <= x^(1/2) < 2^16}
      sum := sum - (Primes16Index(z div Primes16[j]) - j + 1);
    end;
  end;
  for i:=c+1 to b do sum := sum - primepi32(x div Primes16[i]);

  primepi32 := lsumphi32(x,a) + sum;
end;


{$ifdef MPC_PurePascal}
{---------------------------------------------------------------------------}
function _spsp32(a,N,d,s: uint32): boolean;
  {-Strong probable prime test for N with base a, N-1=2^s*d, s>0}
var
  dk,p,N1: uint32;
  i: word;
begin
  _spsp32 := true;
  N1 := N-1;
  {d>0 is odd, so first iteration can be done manually without MulMod}
  p := a;
  dk := d shr 1;
  while dk<>0 do begin
    a := uint64(a)*a mod N;
    if odd(dk) then p := uint64(p)*a mod N;
    dk := dk shr 1;
  end;
  {here p=bases[k]^d mod N}
  if (p=1) or (p=N1) then exit;
  {calculate p^i for i<s, note shifted for loop variable}
  for i:=2 to s do begin
    p := uint64(p)*p mod N;
    if p=N1 then exit;
  end;
  _spsp32 := false;
end;


{-----------------------------------------------------------------------------}
function is_spsp32A(N: longint; const bases: array of longint): boolean;
  {-Strong probable prime test for N with a number of bases. }
  { Return true if N is a probable prime for all bases, or if}
  { N=2. Negative numbers are treated as unsigned (=cardinal)}
var
  s,d,N1: uint32;
  k: integer;
begin
  N1 := uint32(N)-1;
  if odd(N1) or (N1=0) then begin
    is_spsp32A := N=2;
    exit;
  end;
  {calculate d*2^s = N-1, this is done once for all bases}
  d := N1;
  s := 0;
  while not odd(d) do begin
    d := d shr 1;
    inc(s);
  end;
  {default no spsp}
  is_spsp32A := false;
  {now loop through all the bases}
  for k:=low(bases) to high(bases) do begin
    {spsp test for bases[k]}
    if not _spsp32(uint32(bases[k]),uint32(N),d,s) then exit;
  end;
  {All tests passed, this is a spsp}
  is_spsp32A := true;
end;

{$else}

{$ifdef BIT32}
{---------------------------------------------------------------------------}
function _spsp32(a,N,d,s: longint): boolean;
  {-Strong probable prime test for N with base a, N-1=2^s*d, s>0}
var
  res: boolean;
begin
  res := true;
  asm
        pushad
        mov   esi,[N]          {esi=N}
        mov   ecx,[a]          {ecx= a mod N}
        cmp   ecx,esi
        jb    @@0
        mov   eax,ecx          {reduce a mod N}
        sub   edx,edx          {to avoid overflow}
        div   esi
        mov   ecx,edx          {ecx=a mod N}
@@0:    mov   ebx,ecx          {ebx=p=a}
        mov   edi,[d]          {edi=dk}
        shr   edi,1            {dk=dk shr 1}
      {while   dk<>0 do }
@@1:    or    edi,edi
        jz    @@3
        mov   eax,ecx          {a=a^2 mod N}
        mul   eax
        div   esi
        mov   ecx,edx
        shr   edi,1            {dk=dk shr 1}
        jnc   @@1
        mov   eax,ebx          {p=p*a mod N}
        mul   ecx
        div   esi
        mov   ebx,edx
        jmp   @@1
      {end while}

@@3:    dec   ebx              {if p=1 then goto done}
        jz    @@5
        inc   ebx
        mov   edi,[N]
        dec   edi              {edi=N1}
        cmp   ebx,edi          {if p=N1 then goto done}
        jz    @@5
        mov   ecx,s            {remember: s > 0}
        dec   ecx              {if s=1 then N not prime}
        jz    @@np
      {for i=2 to s do}
@@4:    mov   eax,ebx          {p=p^2 mod N}
        mul   ebx
        div   esi
        mov   ebx,edx
        cmp   ebx,edi          {if p=N1 then goto done}
        jz    @@5
        dec   ecx
        jnz   @@4
      {end for}
@@np:   mov   [res],cl         {not prime, here cl=0!}
@@5:
        popad
  end;
  _spsp32 := res;
end;
{$else}
{---------------------------------------------------------------------------}
function _spsp32(a,N,d,s: longint): boolean;
  {-Strong probable prime test for N with base a, N-1=2^s*d, s>0}
var
  res: boolean;
begin
  res := true;
  asm
        db $66;  mov  si,word ptr [N]    {si=N}
        db $66;  mov  cx,word ptr [a]    {cx=a mod N}
        db $66;  cmp  cx,si
                 jb   @@0
        db $66;  mov  ax,cx
        db $66;  sub  dx,dx              {to avoid overflow}
        db $66;  div  si
        db $66;  mov  cx,dx              {cx=a mod N fixed in V1.2.12}
@@0:    db $66;  mov  bx,cx              {bx=p=a}
        db $66;  mov  di,word ptr [d]    {di=dk}
        db $66;  shr  di,1               {dk=dk shr 1}
        {while   dk<>0 do }
@@1:    db $66;  or   di,di
                 jz   @@3
        db $66;  mov  ax,cx              {a=a^2 mod N}
        db $66;  mul  ax
        db $66;  div  si
        db $66;  mov  cx,dx
        db $66;  shr  di,1               {dk=dk shr 1}
                 jnc  @@1
        db $66;  mov  ax,bx              {p=p*a mod N}
        db $66;  mul  cx
        db $66;  div  si
        db $66;  mov  bx,dx
                 jmp  @@1
        {end while}

@@3:    db $66;  dec  bx                 {if p=1 then goto done}
                 jz   @@5
        db $66;  inc  bx
        db $66;  mov  di,word ptr [N]    {di=N1}
        db $66;  dec  di
        db $66;  cmp  bx,di              {if p=N1 then goto done}
                 jz   @@5
                 mov  cx,word ptr [s]    {remember: s > 0}
                 dec  cx                 {if s=1 then N not prime}
                 jz   @@np
        {for i=2 to s do}
@@4:    db $66;  mov  ax,bx              {p=p^2 mod N}
        db $66;  mul  bx
        db $66;  div  si
        db $66;  mov  bx,dx
        db $66;  cmp  bx,di              {if p=N1 then goto done}
                 jz   @@5
                 dec  cx
                 jnz  @@4
        {end for}
@@np:   mov      [res],cl                {not prime, here cl=0!}
@@5:
  end;
  _spsp32 := res;
end;

{$endif}

{-----------------------------------------------------------------------------}
function is_spsp32A(N: longint; const bases: array of longint): boolean;
  {-Strong probable prime test for N with a number of bases. }
  { Return true if N is a probable prime for all bases, or if}
  { N=2. Negative numbers are treated as unsigned (=cardinal)}
var
  s,d,N1: longint;
  k: integer;
begin
  N1 := N-1;
  if odd(N1) or (N1=0) then begin
    is_spsp32A := N=2;
    exit;
  end;
  {calculate d*2^s = N-1, this is done once for all bases}
  {$ifdef BIT16}
    asm
      db $66;  mov  ax,word ptr [N1]
      db $66,       $0F,$BC,$C8     {bsf  ecx,eax}
      db $66;  mov  word ptr [s],cx
      db $66;  shr  ax,cl
      db $66;  mov  word ptr [d],ax
    end;
  {$else}
    asm
      mov  eax,[N1]
      bsf  ecx,eax
      mov  [s],ecx
      shr  eax,cl
      mov  [d],eax
    end;
  {$endif}
  {default no spsp}
  is_spsp32A := false;
  {now loop through all the bases}
  for k:=low(bases) to high(bases) do begin
    {spsp test for bases[k]}
    if not _spsp32(bases[k],N,d,s) then exit;
  end;
  {All tests passed, this is a spsp}
  is_spsp32A := true;
end;
{$endif}


{---------------------------------------------------------------------------}
function IsPrime32(N: longint): boolean;
  {-Test if longint (DWORD) N is prime}
const
  a1: array[0..1] of longint = (31,73);
  a2: array[0..2] of longint = (2,7,61);
type
  LH = packed record L,H: word; end;
begin
  if LH(N).H=0 then begin
    {$ifndef BIT16}
      IsPrime32 := (N=2) or (pbits16[N shr 4] and pmask16[N and $0F] <> 0);
    {$else}
      IsPrime32 := (word(N)=2) or (pbits16[word(N) shr 4] and pmask16[word(N) and $0F] <> 0);
    {$endif}
  end
  else begin
    {Normally here N>=2^16, IN ANY CASE: N MUST NOT BE IN [2,7,31,61,73]!}
    if odd(N) then begin
      {First test N and $80000000 <> 0 selects N>MaxLongint,}
      {second test is the upper bound for spsp set (31,73)] }
      if (N and $80000000 <> 0) or (N>=9080191) then IsPrime32 := is_spsp32A(N, a2)
      else IsPrime32 := is_spsp32A(N, a1);
    end
    else begin
      {N is even, N<>2 therefore N not prime}
      IsPrime32 := false;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function is_primepower32(n: longint; var p: longint): integer;
  {-Test if n is a prime power: return 0 if not, n = p^result otherwise.}
  { Note: contrary to mp_is_primepower a prime n is returned as n = n^1!}
var
  i,r: integer;
begin
  is_primepower32 := 0;
  if n>1 then begin
    if n and 1 = 0 then begin
      {even n can only be powers of 2}
      if n and (n-1) = 0 then begin
        p := 2;
        is_primepower32 := bitsize32(n)-1;
      end;
    end
    else if isprime32(n) then begin
      p := n;
      is_primepower32 := 1;
    end
    else begin
      for i:=2 to NumPrimes16 do begin
        p := Primes16[i];
        if n mod p = 0 then begin
          r := 1;
          n := n div p;
          while (n >= p) and (n mod p = 0) do begin
            n := n div p;
            inc(r);
          end;
          if n=1 then is_primepower32 := r;
          exit;
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function is_spsp32(N, A: longint): boolean;
  {-Strong probable prime test for N with base A, calls is_spsp32A}
begin
  is_spsp32 := is_spsp32A(N, A);
end;

{$ifdef MPC_PurePascal}
  {------------------------------------------------------------}
  function modnpd2(n: longint): integer; {$ifdef HAS_INLINE} inline;{$endif}
    {-Calculate (n mod NPRC_MOD) div 2 (DWORD interpretation)}
  begin
    modnpd2 := (uint32(n) mod NPRC_MOD) shr 1;
  end;
{$else}
  {$ifdef BIT32}
  {------------------------------------------------------------}
  function modnpd2(n: longint): integer; assembler; {&frame-}
    {-Calculate (n mod NPRC_MOD) div 2 (DWORD interpretation)}
  asm
    {$ifdef LoadArgs}
      mov eax,[n]
    {$endif}
     mov  ecx,NPRC_MOD
     sub  edx,edx
     div  ecx
     mov  eax,edx
     shr  eax,1
  end;
  {$else}
  {------------------------------------------------------------}
  function modnpd2(n: longint): integer; assembler;
    {-Calculate (n mod NPRC_MOD) div 2 (DWORD interpretation)}
  asm
    db $66; mov ax, word ptr [n]
    db $66; sub dx,dx
    db $66; mov cx,NPRC_MOD; dw 0;
    db $66; div cx
            mov ax,dx
            shr ax,1
  end;
  {$endif}
{$endif}


{---------------------------------------------------------------------------}
function nextprime32(n: longint): longint;
  {-Next 32 bit prime >= n (DWORD interpretation, see note)}
var
  id,k: integer;
begin
  {easy outs and assure valid range for mod MP_MOD calculation}

  {note: (n>=-4) CANNOT be omitted because of DWORD interpretation}
  {for m=-4=$FFFFFFFC=4294967292, nextprime(m) is greater 2^32 and}
  {will be set to 0 (nextprime(m) mod 2^32 = 15 would be strange!)}

  if (n>=-4) and (n<=7) then begin
    if n<0 then nextprime32 := 0
    else if n<=2 then nextprime32 := 2
    else if n<=3 then nextprime32 := 3
    else if n<=5 then nextprime32 := 5
    else  nextprime32 := 7;
    exit;
  end;

  {make n odd}
  if n and 1 = 0 then inc(n);

  {$ifndef BIT16}
    {avoid warning, id WILL always be initialized (for bug-free modnpd2!)}
    id := 0;
  {$endif}

  {move n to next prime residue class mod MP_MOD and index id into diff array}
  for k:=modnpd2(n) to (NPRC_MOD div 2)-1 do begin
    id := NPRC_OddIdx[k];
    {note: loop terminates via break because NPRC_OddIdx[(NPRC_MOD div 2)-1]<>-1}
    if id<>-1 then break;
    inc(n,2);
  end;

  repeat
    {loop through possible primes}
    if IsPrime32(n) then begin
      nextprime32 := n;
      exit;
    end;
    {move to next candidate}
    inc(n,longint(NPRC_Diff[id]));
    {get next increment index}
    inc(id); if id>=NPRC_NRC then id:=0;
  until false;
end;


{---------------------------------------------------------------------------}
procedure FindFirstPrime32(n: longint; var ctx: TPrimeContext);
  {-Find first prime >= n and initialize ctx for FindNextPrime32}
begin
  with ctx do begin
    next  := n;
    index := -1;
  end;
  FindNextPrime32(ctx);
end;


{---------------------------------------------------------------------------}
procedure FindNextPrime32(var ctx: TPrimeContext);
  {-Find next 32 bit prime (DWORD interpretation, see note), success if ctx.prime<>0}
var
  k: integer;
  found: boolean;
const
  MP32 = longint($FFFFFFFB); {4294967291 = prevprime(2^32)}
begin
  with ctx do begin
    if index<0 then begin
      {note: (n>=-4) CANNOT be omitted because of DWORD interpretation}
      if (next>=-4) and (next<=7) then begin
        if next<0 then begin
          prime := 0;
          next  := -1;
        end
        else if next<=2 then begin
          prime := 2;
          next  := 3;
        end
        else if next<=3 then begin
          prime := 3;
          next  := 5;
        end
        else if next<=5 then begin
          prime := 5;
          next  := 7;
        end
        else  begin
          prime := 7;
          next  := 11;
        end;
        exit;
      end;
      {first index calculation after FindFirstPrim32}
      {make n odd}
      if next and 1 = 0 then inc(next);
      {move n to next prime residue class mod MP_MOD and index id into diff array}
      for k:=modnpd2(next) to (NPRC_MOD div 2)-1 do begin
        index := NPRC_OddIdx[k];
        {note: loop terminates via break because NPRC_OddIdx[(NPRC_MOD div 2)-1]<>-1}
        if index<>-1 then break;
        inc(next,2);
      end;
    end;
    repeat
      {loop through possible primes}
      found := IsPrime32(next);
      if found then begin
        prime := next;
        if next=MP32 then begin
          next  := -1;
          index := -2;
          exit;
        end;
      end;
      {move to next candidate}
      inc(next,longint(NPRC_Diff[index]));
      {get next increment index}
      inc(index);
      if index>=NPRC_NRC then index:=0;
    until found;
  end;
end;


{---------------------------------------------------------------------------}
procedure nextprime32_array(n,cmax: longint; var a: array of longint);
  {-Fill an array with the next 32 bit primes >= n (DWORD interpretation).}
  { If a[i0]=4294967291 then a[i]=0 for i>i0. If cmax<=0 fill complete array.}
var
  i,k,ma: longint;
  ctx: TPrimeContext;
begin
  ma := high(a);
  if cmax>0 then begin
    dec(cmax);
    if cmax<ma then ma := cmax;
  end;
  FindFirstPrime32(n, ctx);
  for i:=low(a) to ma do begin
    a[i] := ctx.prime;
    if ctx.prime=0 then begin
      {if no more 32 bit primes fill rest of array with 0}
      for k:=i+1 to ma do a[k] := 0;
      break;
    end;
    FindNextPrime32(ctx);
  end;
end;


{---------------------------------------------------------------------------}
function prevprime32(n: longint): longint;
  {-Previous 32 bit prime <= n, prevprime32(0)=0, (DWORD interpretation)}
var
  id,k: integer;
begin
  {easy outs and assure valid range for mod MP_MOD calculation}
  {note: (n>=0) CANNOT be omitted because of DWORD interpretation}

  if (n>=0) and (n<11) then begin
    if n<2 then prevprime32 := 0
    else if n<3 then prevprime32 := 2
    else if n<5 then prevprime32 := 3
    else if n<7 then prevprime32 := 5
    else prevprime32 := 7;
    exit;
  end;

  {make n odd}
  if n and 1 = 0 then dec(n);

  {$ifndef BIT16}
    {avoid warning, id WILL always be initialized (for bug-free modnpd2!)}
    id := 0;
  {$endif}

  {move n to prev prime residue class mod MP_MOD and index id into diff array}
  for k:=modnpd2(n) downto 0 do begin
    id := NPRC_OddIdx[k];
    {note: loop is always terminated via break because NPRC_OddIdx[0]<>-1}
    if id<>-1 then break;
    dec(n,2);
  end;

  repeat
    {loop through possible primes}
    if IsPrime32(n) then begin
      prevprime32 := n;
      exit;
    end;
    {get prev increment index}
    dec(id); if id<0 then id:=NPRC_NRC-1;
    {move to prev candidate}
    dec(n,longint(NPRC_Diff[id]));
  until false;
end;


{---------------------------------------------------------------------------}
function safeprime32(n: longint): longint;
  {-Return the next safe prime p >= n, i.e. p and (p-1)/2 prime; 0 if n > 2147483579}
var
  ctx: TPrimeContext;
  p: longint;
begin
  safeprime32 := 5;
  if n>5 then begin
    if n>=2147483579 then begin
      if n>2147483579 then safeprime32 := 0
      else safeprime32 := n;
      exit;
    end;
    FindFirstPrime32(n, ctx);
    repeat
      p := (ctx.prime-1) shr 1;
      if odd(p) and IsPrime32(p) then begin
        safeprime32 := ctx.prime;
        exit;
      end;
      FindNextPrime32(ctx);
    until false;
  end;
end;

{---------------------------------------------------------------------------}
{--------------------- 32-bit number-theoretic functions -------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
procedure PrimeFactor32(n: longint; var FR: TFactors32);
  {-Return the prime factorization of n > 1, FR.pcount=0 if n < 2}
var
  i: integer;
  p: word;
  e: byte;
type
  LH = packed record L,H: word; end;
begin
  fillchar(FR, sizeof(FR),0);
  with FR do begin
    FR.number := n;
    if n<2 then exit;
    {Get powers of two}
    if integer(n) and 1 = 0 then begin
      e := 0;
      pcount := 1;
      primes[1] := 2;
      repeat
        inc(e);
        n := n shr 1;
      until integer(n) and 1 <> 0;
      pexpo[1] := e;
      if n=1 then exit;
    end;
    {Fast check if remainder is small prime}
    if (LH(n).H=0) and IsPrime16(LH(n).L) then begin
      inc(pcount);
      pexpo[pcount]  := 1;
      primes[pcount] := n;
      exit;
    end;
    {Process primes up to sqrt(remaining part)}
    for i:=2 to Primes16Index(isqrt32(n)) do begin
      p := Primes16[i];
      if n mod p = 0 then begin
        e := 1;
        n := n div p;
        while n mod p = 0 do begin
          inc(e);
          n := n div p;
        end;
        inc(pcount);
        {$ifdef MPC_USE_Assert}
          {only if assert supported by compiler or debug}
          assert(pcount<=10, MPAF+'pcount<=10 in PrimeFactor32');
        {$endif}
        pexpo[pcount]  := e;
        primes[pcount] := p;
        if n=1 then exit;
      end;
    end;
    if n>1 then begin
      {Remaining factor of n must be prime}
      inc(pcount);
      {$ifdef MPC_USE_Assert}
        {only if assert supported by compiler or debug}
        assert(pcount<=10, MPAF+'pcount<=10 in PrimeFactor32');
      {$endif}
      pexpo[pcount]  := 1;
      primes[pcount] := n;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function EulerPhi32(n: longint): longint;
  {-Return Euler's totient function phi(|n|), phi(0)=0. For n > 0 }
  { this the number of positive integers k <= n with gcd(n,k) = 1.}
var
  FR: TFactors32;
  i,j: integer;
  f,p: longint;
begin
  n := abs(n);
  if n<2 then EulerPhi32 := n
  else with FR do begin
    PrimeFactor32(n,FR);
    f := 1;
    {Just the product: (primes[i]-1)*primes[i]^(pexpo[i]-1)}
    for i:=1 to pcount do begin
      p := primes[i];
      f := f*(p-1);
      for j:=2 to pexpo[i] do f := f*p;
    end;
    EulerPhi32 := f;
  end;
end;


{---------------------------------------------------------------------------}
function Carmichael32(n: longint): longint;
  {-Return the Carmichael function lambda(|n|), lambda(0)=0. For n > 0 this}
  { is the least k such that x^k = 1 (mod n) for all x with gcd(x,n) = 1.  }
var
  FR: TFactors32;
  i,j,k: integer;
  f,p,c: longint;
const
  CS: array[0..4] of byte = (0,1,1,2,2);
begin
  n := abs(n);
  if n<=4 then Carmichael32 := CS[n]
  else with FR do begin
    PrimeFactor32(n,FR);
    c := 1;
    for i:=1 to pcount do begin
      p := primes[i];
      {compute f = lambda(primes[i]^pexpo[i])}
      f := p-1;
      k := pexpo[i];
      if k>1 then begin
        for j:=2 to k do f := f*p;
        if (p=2) and (k>2) then f := f shr 1;
      end;
      if c=1 then c := f
      else begin
        {c := lcm(c,f)}
        p := gcd32(f,c);
        if p>1 then c := c div p;
        c := f*c;
      end;
    end;
    Carmichael32 := c;
  end;
end;


{---------------------------------------------------------------------------}
function is_Carmichael32(n: longint): boolean;
  {-Test if |n| is a Carmichael number}
var
  i,cnt: integer;
  p: word;
  nm1: longint;
type
  LH = packed record L,H: word; end;
begin
  {Checked against https://oeis.org/A002997/b002997.txt for all n}
  n := abs(n);
  is_Carmichael32 := false;
  {Facts:
    a) [10, Theorem 3.4.6 (Korselt criterion)]: An integer n is a Carmichael
       number if and only if n is positive, composite, squarefree, and for
       each prime p dividing n we have p-1 dividing n-1.
    b) [7, Ch.2 IX] Every Carmichael number is odd and is the product of
       three or more distinct primes. The smallest Carmichael number is 561.
  }
  {Exclude small or even n}
  if (n and 1 = 0) or (n<561) then exit;

  {Fast check if remainder is small prime}
  if (LH(n).H=0) and IsPrime16(LH(n).L) then exit;

  {Process primes up to sqrt(remaining part)}
  cnt := 0;
  nm1 := n-1;
  for i:=2 to Primes16Index(isqrt32(n)) do begin
    p := Primes16[i];
    if n mod p = 0 then begin
      {check p-1|n-1}
      if nm1 mod (p-1) <> 0 then exit;
      {check squarefree}
      n := n div p;
      if n mod p = 0 then exit;
      inc(cnt);
      if n=1 then break;
    end;
  end;
  if n>1 then begin
    {Remaining factor of n must be prime}
    if (nm1 mod (n-1) <> 0) then exit;
    inc(cnt);
  end;
  if cnt > 2 then is_Carmichael32 := true;
end;


{---------------------------------------------------------------------------}
function Moebius32(n: longint): integer;
  {-Return the Moebius function mu(abs(n)), mu(0)=0, mu(1)=1. mu(n)=(-1)^k }
  { if n > 1 is the product of k different primes; otherwise mu(n)=0.      }
var
  i,k,m: integer;
  q: longint;
  p: word;
const
  m07: array[0..7] of integer = (0,1,-1,-1,0,-1,1,-1);
type
  LH = packed record L,H: word; end;
begin
  n := abs(n);
  if n<8 then begin
    Moebius32 := m07[n];
    exit;
  end;
  {Set default value}
  Moebius32 := 0;
  case n and 3 of
      0: exit; {multiple of 4}
      2: begin
           n := n shr 1;
           if n<8 then begin
             Moebius32 := -m07[n];
             exit;
           end;
           k := 1;
         end;
    else k := 0;
  end;

  {k counts the number of different prime factors}
  if (LH(n).H=0) and IsPrime16(LH(n).L) then begin
    {remaining part is a 16-bit prime}
    inc(k);
  end
  else begin
    {Let n0 current n; check primes up to n0^(1/3)}
    m := icbrt32(n); {n>=8, i.e. m>=2}
    for i:=2 to NumPrimes16 do begin
      p := Primes16[i];
      if p>m then break;
      q := n div p;
      if q<p then break; {done and no factor/square}
      if p*q=n then begin
        {p is factor}
        if q mod p = 0 then exit; {p^2 | n}
        {continue with p removed from n}
        n := q;
        inc(k);
      end;
    end;
    if n>1 then begin
      {The prime factors of the remaining current n are greater than  }
      {n0^(1/3), i.e. there are at most two remaining factors: n is   }
      {either a square, a prime, or a product of two different primes.}
      if isprime32(n) then inc(k)
      else if is_square32(n) then exit
      else inc(k,2);
    end;
  end;
  if odd(k) then Moebius32 := -1
  else Moebius32 := 1;
end;


{---------------------------------------------------------------------------}
function is_primroot32(const g,n: longint): boolean;
  {-Test if g is primitive root mod n}
var
  p,phi: longint;
  i,j: integer;
  FR: TFactors32;
begin
  {Handle trivial and small cases}
  is_primroot32 := false;
  if n<5 then begin
    if n>1 then is_primroot32 := g=n-1;
    exit;
  end;
  if (n and 3 = 0) or (is_square32(g)) or (gcd32(g,n)<>1) then exit;

  {Here n>4, n mod 4 <> 0. Use prime factorizations of n and phi(n)}
  PrimeFactor32(n,FR);
  with FR do begin
    {exit if at least two odd prime factors}
    if (pcount>2) or ((pcount=2) and odd(primes[1]) and odd(primes[2])) then exit;
    {Compute phi(n). Here n = 2*p^k or p^k, so phi(n) = phi(p^k)}
    i := 1;
    if (pcount>1) and (primes[1]=2) then inc(i);
    p := primes[i];
    phi := p-1;
    for j:=2 to pexpo[i] do phi := phi*p;
  end;

  PrimeFactor32(phi,FR);
  with FR do begin
    {Check if g is a primitive root using Adler/Coury[26], Theorem 6.8:}
    {"If (g,n) = 1, then g is a primitive root of n if and only if     }
    { g^(phi(n)/q) <> 1 mod n for every prime divisor q of phi(n)."    }
    for i:=1 to pcount do begin
      if exptmod32(g, phi div primes[i], n) = 1 then exit;
    end;
  end;
  is_primroot32 := true;
end;


{---------------------------------------------------------------------------}
function is_squarefree32(n: longint): boolean;
  {-Return true if n is squarefree, i.e. not divisible by a square > 1}
type
  LH = packed record L,H: word; end;
var
  i: integer;
  q: longint;
  p,m: word;
begin
  n := abs(n);
  if n<4 then begin
    is_squarefree32 := n<>0;
    exit;
  end;
  {Set default value}
  is_squarefree32 := false;
  case n and 3 of
    0: exit;            {multiple of 4}
    2: n := n shr 1;    {single factor 2, continue with n/2}
  end;

  if (LH(n).H=0) and IsPrime16(LH(n).L) then begin
    {remaining part is a 16-bit prime}
    is_squarefree32 := true;
    exit;
  end;

  {Check primes up to n^(1/3)}
  m := icbrt32(n);
  for i:=2 to NumPrimes16 do begin
    p := Primes16[i];
    if p>m then break;
    q := n div p;
    if q<p then begin
      is_squarefree32 := true;
      exit;
    end;
    if p*q=n then begin
      {p is factor}
      if q mod p = 0 then exit;
      {continue with p removed from n}
      n := q;
    end;
  end;
  {Primes <= n^(1/3) are removed, so current n is prime or has two factors}
  is_squarefree32 := not is_square32(n);
end;


{---------------------------------------------------------------------------}
function is_fundamental32(d: longint): boolean;
  {-Return true, if d is a fundamental discriminant (either d=1 mod 4 and }
  { squarefree, or d=0 mod 4, d/4 = 2,3 mod and squarefree), false if not.}
var
  r: integer;
  a: longint;
begin
  {Ref: Pari[12], function isfundamental in arith1.c}
  is_fundamental32 := false;
  a := abs(d);
  r := a and 15;
  if r=0 then exit;
  if r and 3 = 0 then begin
    r := r shr 2;
    if d<0 then r := 4-r;
    if r<>1 then is_fundamental32 := is_squarefree32(a shr 2);
  end
  else begin
    r := r and 3;
    if d<0 then r := 4-r;
    if r=1 then is_fundamental32 := is_squarefree32(a);
  end;
end;


{---------------------------------------------------------------------------}
function order32(n,m: longint): longint;
  {-Return the order of n mod m, m > 1, i.e. the smallest integer}
  { e with n^e = 1 mod m; if gcd(n,m) <> 1 the result is 0.}
var
  i,j,k: integer;
  e,p,y: longint;
  FR: TFactors32;
begin
  if (m<2) or (gcd32(n,m)<>1) then begin
    order32 := 0;
    exit;
  end;
  if (n<0) or (n>=m) then begin
    n := n mod m;
    {Fix possible Pascal 'feature'}
    if n<0 then n := n+m;
  end;

  e := Carmichael32(m);
  PrimeFactor32(e,FR);

  {Algorithm 1.4.3 from Cohen[24], see also HAC[5], 4.79}
  with FR do begin
    for i:=1 to pcount do begin
      k := pexpo[i];
      if k>1 then begin
        p := primes[i];
        {compute e := e / primes[i] ^ pexpo[i] and y = n^e mod m}
        if p=2 then e := e shr k
        else begin
          y := p;
          for j:=pred(k) downto 1 do y := y*p;
          e := e div y;
        end;
        y := exptmod32(n,e,m);
        {compute local order}
        while y<>1 do begin
          e := e*p;
          y := exptmod32(y, p, m);
        end;
      end
      else begin
        {k=1: same as above with full simplification}
        y := e div primes[i];
        if exptmod32(n,y,m)=1 then e := y;
      end
    end;
    order32 := e;
  end;
end;


{---------------------------------------------------------------------------}
function core32(n: longint): longint;
  {-Return the squarefree part of n, n<>0, i.e. the unique squarefree integer c with n=c*f^2}
var
  c: longint;
  i: integer;
  FR: TFactors32;
begin
  {define core32(0)=0}
  if n=0 then c:=0
  else begin
    PrimeFactor32(abs(n),FR);
    with FR do begin
      if n<0 then c := -1 else c := 1;
      for i:=1 to pcount do begin
        if odd(pexpo[i]) then c := c*primes[i];
      end;
    end;
  end;
  core32 := c;
end;


{---------------------------------------------------------------------------}
function rad32(n: longint): longint;
  {-Return the radical rad(n) of n = product of the distinct prime factors of n.}
var
  r: longint;
  i: integer;
  FR: TFactors32;
begin
  {define rad32(0)=0}
  if n=0 then r:=0
  else begin
    PrimeFactor32(abs(n),FR);
    with FR do begin
      if n<0 then r := -1 else r := 1;
      for i:=1 to pcount do r := r*primes[i];
    end;
  end;
  rad32 := r;
end;


{---------------------------------------------------------------------------}
function quaddisc32(n: longint): longint;
  {-Return the discriminant of the quadratic field Q(sqrt(n))}
var
  c: longint;
  r: integer;
begin
  {Ref: Pari[12], function quaddisc in arith1.c}
  c := core32(n);
  r := abs(c) and 3;
  if n<0 then r := 4-r;
  if r>1 then c := 4*c;
  quaddisc32 := c;
end;


{---------------------------------------------------------------------------}
function prime32(k: longint): longint;
  {-Return the kth prime if 1 <= k <= 105097565, 0 otherwise}
var
  i,pk: longint;
  c,l1,l2: double;
  ctx: TPrimeContext;
const
  c1 = 2.006218;
  c2 = 2.004383;
  k1 = 57000000;
begin
  if (k<1) or (k>105097565) then prime32 := 0
  else if k<=6542 then prime32 := Primes16[k]
  else begin
    {This is the slightly modified asymptotic expansion for pk (where c=2), see}
    {e.g. http://functions.wolfram.com/13.03.06.0004.01 or P. Dusart, 1999, The}
    {kth prime is greater than k(log k + log log k - 1) for k>=2, Math.Comp.68,}
    {p.411. http://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-01037-6/}
    if k<k1 then c := c1 else c := c2;
    l1 := ln(k);
    l2 := ln(l1);
    {Get approximation for pk, compute primepi, and correct}
    pk := trunc(k*(l1+l2-1.0 + (l2-c)/l1 - ((l2-6.0)*l2 + 11.0)/(2.0*l1*l1)));
    {pk <= MaxLongint; verified by exhaustive computation with the procedure }
    {check_pk_formula4. The maximum deviation of 42198 occurs for k=95550393,}
    {this gives an i value from primepi32 which is 2314 to low.}
    FindFirstPrime32(pk,ctx);
    i := primepi32(ctx.prime);
    while i<k do begin
      FindNextPrime32(ctx);
      inc(i);
    end;
    prime32 := ctx.prime;
  end;
end;


{---------------------------------------------------------------------------}
function primroot32(n: longint): longint;
  {-Compute the smallest primitive root mod n, 0 if n does not have a prim.root}
var
  g,phi: longint;
  i,j: integer;
  FR: TFactors32;
  nprim, neven: boolean;
const
  spr: array[2..31] of byte = (1,2,3,2,5,3,0,2,3,2,0,2,3,0,0,3,5,2,0,0,7,5,0,2,7,2,0,2,0,3);
begin
  {Handle trivial and small cases}
  primroot32 := 0;
  if n<2 then exit;
  if n<32 then begin
    primroot32 := spr[n];
    exit;
  end;
  {Check multiple of 4 AFTER the small cases, otherwise primroot32(4) would be 0}
  if n and 3 = 0 then exit;

  {here n>32, use prime factorizations of n and phi(n)}
  PrimeFactor32(n,FR);
  with FR do begin
    {exit if at least two odd prime factors}
    if (pcount>2) or ((pcount=2) and odd(primes[1]) and odd(primes[2])) then exit;
    {Remember if n=p is prime: Avoid gcd calculations in search loop}
    nprim := (pcount=1) and (pexpo[1]=1);
    neven := not odd(n);
    {Compute phi(n). Here n = 2*p^k or p^k, so phi(n) = phi(p^k)}
    i := 1;
    if (pcount>1) and (primes[1]=2) then inc(i);
    g := primes[i];
    phi := g-1;
    for j:=2 to pexpo[i] do phi := phi*g;
  end;

  PrimeFactor32(phi,FR);
  for g:=2 to n-1 do begin
    {quick check gcd(g,n)=1}
    if neven and (g and 1 = 0) then continue;
    {exclude squares which are never primitive roots}
    if is_square32(g) then continue;
    if nprim or (gcd32(g,n)=1) then with FR do begin
      {Check if g is a primitive root using Adler/Coury[26], Theorem 6.8:}
      {"If (g,n) = 1, then g is a primitive root of n if and only if     }
      { g^(phi(n)/q) <> 1 mod n for every prime divisor q of phi(n)."    }
      i := 1;
      while i <= pcount do begin
        if exptmod32(g, phi div primes[i], n) = 1 then break;
        inc(i);
      end;
      if i>pcount then begin
        primroot32 := g;
        exit;
      end;
    end;
  end;
end;

(*
{---------------------------------------------------------------------------}
function primroot32p(p: longint): longint;
  {-Compute the smallest primitive root mod p, p odd prime}
var
  g,h,p1: longint;
  i,k: integer;
  FR: TFactors32;
begin

  primroot32p := 0;
  p1 := p-1;
  if odd(p1) then exit;

  PrimeFactor32(p1,FR);
  h := 1;
  for g:=2 to p1 do begin
    {exclude squares which are never primitive roots}
    if is_square32(g) then continue;
    if exptmod32(g, p1, p) <> 1 then exit;
    with FR do begin
      for i:=1 to pcount do begin
        h := exptmod32(g, p1 div primes[i], p);
        if h=1 then break;
      end;
    end;
    if h<>1 then begin
      primroot32p := g;
      exit;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function primroot32fast(n: longint): longint;
  {-Compute a primitive root mod n, 0 if n does not have a primitive root. }
  { If n is not prime, the result is not necessarily the smallest primroot.}
var
  g,p,p2,pk: longint;
  FR: TFactors32;
begin
  if (n=2) or (n=4) then begin
    primroot32fast := n-1;
    exit;
  end;
  primroot32fast := 0;
  if (n<2) or (n and 3 = 0) then exit;
  if odd(n) then pk := n
  else pk := n shr 1;
  {here p is odd and > 2}
  PrimeFactor32(pk,FR);
  if FR.pcount<>1 then exit; {product of at least two odd primes}
  p := FR.primes[1];
  g := primroot32p(p);
  if g=0 then exit;  {should not happen}
  if FR.pexpo[1] > 1 then begin
    {case p^k, g is primitive root mod p. Check if g is a primroot mod p^2.}
    {If yes then g is a primroot mod p^k, if not, then g + p is a primroot }
    {mod p^k. See Cohen[24], Lemma 1.4.5 and Forster[4], Satz 8.11.        }
    p2 := sqr(p);
    if exptmod32(g,p-1,p2) <> 1 then g := g + p mod p2;
  end;
  if (not odd(n)) and (not odd(g)) then begin
    {Case n = 2*p^k: See Cohen[24] last sentence in section 1.4.1}
    g := g + pk;
  end;
  primroot32fast := g;
end;
*)

{---------------------------------------------------------------------------}
function tau32(n: longint): integer;
  {-Return the number of positive divisors of n (including 1 and n)}
var
  i,t: integer;  {t_max=1600 for n=2095133040}
  FR: TFactors32;
begin
  if n=0 then tau32 := 0
  else begin
    PrimeFactor32(abs(n), FR);
    with FR do begin
      t := 1;
      for i:=1 to pcount do t := t*(1+pexpo[i]);
    end;
    tau32 := t;
  end;
end;


{---------------------------------------------------------------------------}
function dlog32_ex(a,b,p,JMAX: longint): longint;
  {-Compute the discrete log_a(b) mod p using Pollard's rho  }
  { algorithm; p prime, a > 1, b > 0; return -1 for failure. }
  { Expanded version with variable trial parameter JMAX >= 0.}
var
  x1,a1,b1,x2,a2,b2,i,p1,p13,p23,q,d: longint;
  j: integer;
const
  IMAX = 215000;  {Max. number of iterations per trial, ~ sqrt(2^31*ln(2^31))}
  MLIM = 46340;   {isqrt(2^31), x*y mod p does not overflow with 32-bit arith}

  {Ref: J.M.Pollard, "Monte Carlo methods for index computation (mod p)",}
  {Mathematics of Computation, 32 (1978), p.918-924. See also R.Crandall,}
  {C.Pomerance [10], Chapter 5.2.2 and HAC [5] Alg. 3.60. Note that these}
  {references assume that a is a primitive root, but here this assumption}
  {is dropped: This is more general but the function may fail.           }

  procedure step(var xi,ai,bi: longint);
    {-Perform next step of pseudo-random sequence}
  begin
    if xi<p13 then begin
      xi := mulmod32(b,xi,p);
      bi := bi+1; if bi>=p1 then dec(bi,p1);
    end
    else if xi<p23 then begin
      xi := mulmod32(xi,xi,p);
      if ai<=q then ai := ai+ai else ai := (ai-p1)+ai;
      if bi<=q then bi := bi+bi else bi := (bi-p1)+bi;
    end
    else begin
      xi := mulmod32(a,xi,p);
      ai := ai+1;  if ai>=p1 then dec(ai,p1);;
    end;
  end;

  procedure linear_search;
    {-linear search, assumes a<>b}
  var
    k: longint;
  begin
    d := a;
    for k:=2 to p1 do begin
      if p < MLIM then d := d*a mod p
      else d := mulmod32(a,d,p);
      if d=b then begin
        dlog32_ex := k;
        exit;
      end
      else if d=a then break;
    end;
    dlog32_ex := -1;
  end;

begin

  {in any case return log_a(1) = 0}
  if b=1 then begin
    dlog32_ex := 0;
    exit;
  end;

  {Default return value}
  dlog32_ex := -1;

  {Normalize and check p,a,b}
  if not isprime32(p) then exit;
  p1 := p-1;
  if a>=p then a := a mod p;
  if b>=p then b := b mod p;
  if b=a then begin
    dlog32_ex := 1;
    exit;
  end;
  if (b<1) or (a<2) or (a=p1) then exit;

  if p < 20 then begin
    {Here we have just to few good collisions and/or much overhead}
    {Note: Benchmarking ALL pairs (a,b), a=2..p-2, b=1..p-2, the  }
    {linear search is faster for safe primes up to 107.           }
    linear_search;
    exit;
  end;

  q := (p-1) div 2;
  if b=p1 then begin
    {b=-1, if a is a primitive root this gives dlog32 = q, so check this.}
    if exptmod32(a,q,p)=b then begin
      dlog32_ex := q;
      exit;
    end;
  end;

  {boundaries for partition}
  p13 := p div 3;
  p23 := 2*p13;

  {Try up to JMAX starting values for cycle finding}
  for j:=0 to JMAX do begin
    if j=0 then begin
      {Standard starting values}
      x1 := 1;
      a1 := 0;
      b1 := 0;
    end
    else if j=1 then begin
      x1 := a;
      a1 := 1;
      b1 := 0;
    end
    else begin
      a1 := 1 + mp_random_long mod (p1-1);
      b1 := 1 + mp_random_long mod (p1-1);
      x1 := exptmod32(a,a1,p);
      x2 := exptmod32(b,b1,p);
      x1 := mulmod32(x1,x2,p);
    end;

    x2 := x1;
    a2 := a1;
    b2 := b1;

    i := 0;
    while i<IMAX do begin
      {Floyd's cycle finding}
      step(x1,a1,b1);
      step(x2,a2,b2);
      step(x2,a2,b2);
      if x1=x2 then break;
      inc(i);
    end;

    if i<IMAX then begin
      {Cycle found: a^a1*b^b1=x1=x2=a^a2*b^b2 mod p}
      {So solve (b1-b2)*dlog_a(b) = (a2-a1) mod p1 }
      {This is only valid if a is a primitive root,}
      {if not, solutions may not exist or be unique}
      b1 := b1-b2;
      a1 := a2-a1;
      if b1<0 then b1 := b1+p1;
      if a1<0 then a1 := a1+p1;
      if b1<>0 then begin
        {With the reduced differences solve b1*x = a1 mod p1}
        xgcd32(b1,p1,x1,i,d);
        if a1 mod d = 0 then begin
          while x1<0 do x1 := x1+p1;
          x1 := mulmod32(x1,a1 div d,p1);
          x2 := p1 div d;
          {solutions are x1+i*x2 mod p1, i=0..d-1}
          for i:=0 to d-1 do begin
            a2 := (x1 + i*x2) mod p1;
            if exptmod32(a,a2,p) = b then begin
              dlog32_ex := a2;
              exit;
            end;
          end;
        end;
      end;
    end;
  end; {end j}

  {No solution found with the randomized algorithm. Do a brute}
  {force linear search if p <= IMAX. With the standard dlog32 }
  {parameters some of the rare cases observed for p > 19 are: }
  {p=179, a=127, b=52; or p=4703, a=3204, b=2318.             }

  if p<=IMAX then begin
    linear_search;
  end;

end;


{---------------------------------------------------------------------------}
function dlog32(a,b,p: longint): longint;
  {-Compute the discrete log_a(b) mod p using Pollard's rho algorithm: i.e.}
  { solve a^x = b mod p, with p prime, a > 1, b > 0; return -1 for failure.}
  { If a is no primitive root mod p, solutions may not exist or be unique. }
const
  IMAX = 215000;  {Max. number of iterations per trial, ~ sqrt(2^31*ln(2^31))}
  MLIM = 46340;   {isqrt(2^31), x*y mod p does not overflow with 32-bit arith}
var
  JMAX: longint;
begin
  {Set number of trials, for small p<IMAX there will be }
  {a final linear search, therefore JMAX can be smaller.}
  if p<MLIM then JMAX := 4
  else if p<IMAX then JMAX := 8
  else JMAX := 32;
  {Note that there are annoying failure cases even for large JMAX, e.g. }
  {b=220285; a=286566; p= 420263. But here a is NO primitive root mod p.}
  {The failure reason is that for most collisions d=2 and a2-a1 is odd! }
  dlog32 := dlog32_ex(a,b,p,JMAX);
end;


{---------------------------------------------------------------------------}
{------------------------ Prime sieve functions ----------------------------}
{---------------------------------------------------------------------------}

{$ifdef MPC_SmallSieve}
  {$ifdef BIT16}
    {$define SieveUseWords}
  {$endif}
{$endif}

{The prime sieve routines are based on Jason Papadopoulos's public domain }
{code prime_sieve.c available from http://www.boo.net/~jasonp/qs.html (in }
{the msieve source archive). Note, that this is a heavily modified version}
{which uses booleans instead of bits. If MPC_SmallSieve is defined, then  }
{the maximum prime is tailored to  1908867043, which is the largest prime }
{below MaxLongint with all paux^[i] offsets <= $ffff; if MPC_SmallSieve is}
{undefined the maximum prime is MaxLongint 2^31-1. Jason generates the 16 }
{bit primes on the fly, numbered from 0. We use the Primes16 array with   }
{Primes16[1]=2, therefore the paux^[i] are indexed with the shift +2. And }
{my procedures do some checks and return 1 if all primes are generated.   }


const
  SIEVE_BLOCK2   = 2*SIEVE_BLOCK;
  SIEVE_MAXBLOCK = 1 + SIEVE_MAXPRIME div SIEVE_BLOCK2;


{-------------------------------------------------------------------------}
procedure prime_sieve_nextblock(var sieve: TSieve);
  {-Perform sieving for the next block}
var
  i: integer;
  {$ifdef SieveUseWords}
    p,r: word;
  {$else}
    p,r: longint;
  {$endif}
begin
  with sieve do begin
    fillchar(psieve^,SIEVE_BLOCK,ord(false));
    for i:=0 to (SIEVE_PRIMES-1) do begin
      p := Primes16[i+2];
      r := paux^[i];
      while r < SIEVE_BLOCK do begin
        psieve^[r] := true;
        inc(r,p);
        {$ifdef SieveUseWords}
          {detect overflow, note that inc/dec do not produce overflow RTE}
          if r<p then break;
        {$endif}
      end;
      dec(r,SIEVE_BLOCK);
      {$ifdef MPC_SmallSieve}
        paux^[i] := word(r);
      {$else}
        paux^[i] := r;
      {$endif}
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure prime_sieve_reset(var sieve: TSieve; first_prime: longint);
  {-Initialize already allocated sieve to return primes >= first_prime; first_prime <= SIEVE_MAXPRIME}
var
  i: integer;
  block_start: longint;
  {$ifdef SieveUseWords}
    p,r: word;
  {$else}
    p,r: longint;
  {$endif}
begin
  with sieve do begin
    if (paux=nil) or (psieve=nil) {$ifdef MPC_SmallSieve}or (first_prime>SIEVE_MAXPRIME){$endif} then begin
      {Out of memory/not allocated or first_prime too large, indicate 'done'}
      curr_blk := SIEVE_MAXBLOCK;
      exit;
    end;

    fillchar(psieve^,SIEVE_BLOCK,ord(false));

    if first_prime and 1 <> 0 then dec(first_prime)
    else if first_prime<=2 then first_prime := 1;

    {Calculate block number and offset of first prime}
    curr_blk :=  first_prime div SIEVE_BLOCK2;
    curr_off := (first_prime mod SIEVE_BLOCK2) shr 1;

    {Calculate the sieve offsets into the current block. The sieve is}
    {compressed so that even multiples of sieving primes are skipped.}
    block_start := curr_blk*SIEVE_BLOCK2;
    if block_start=0 then begin
      {if preparing block 0, also skip the first sieve update}
      for i:= 0 to (SIEVE_PRIMES-1) do begin
        p := Primes16[i+2];
        paux^[i] := p + p shr 1;
      end;
    end
    else begin
      for i:=0 to (SIEVE_PRIMES-1) do begin
        p := Primes16[i+2];
        r := p - (block_start mod p);
        if r and 1 = 0 then inc(r,p);
        paux^[i] := r shr 1;
      end;
    end;
  end;
  prime_sieve_nextblock(sieve);
end;


{---------------------------------------------------------------------------}
procedure prime_sieve_init(var sieve: TSieve; first_prime: longint);
  {-Allocate/initialize sieve to return primes >= first_prime; first_prime <= SIEVE_MAXPRIME-SIEVE_BLOCK}
begin
  with sieve do begin
    paux   := mp_getmem(sizeof(TAux));
    psieve := mp_getmem(sizeof(TFlags));
    if (paux=nil) or (psieve=nil) then begin
      {Out of memory, indicate 'done'}
      curr_blk := SIEVE_MAXBLOCK;
      exit;
    end;
  end;
  prime_sieve_reset(sieve, first_prime);
end;


{---------------------------------------------------------------------------}
procedure prime_sieve_clear(var sieve: TSieve);
  {-Release memory allocated by prime_sieve_init}
begin
  with sieve do begin
    mp_freemem(pointer(paux),sizeof(TAux));
    mp_freemem(pointer(psieve),sizeof(TFlags));
  end;
end;


{---------------------------------------------------------------------------}
function prime_sieve_next(var sieve: TSieve): longint;
  {-Return next prime from sieve, 1 if done}
var
  {$ifdef SieveUseWords}
    off: word;
  {$else}
    off: longint;
  {$endif}
begin
  with sieve do begin
    {$ifdef debug}
      assert((psieve<>nil) and (paux<>nil), MPAF+'(psieve<>nil) and (paux<>nil)');
    {$endif}
    if curr_blk>=SIEVE_MAXBLOCK then begin
      prime_sieve_next := 1;
      exit;
    end;
    off := curr_off;
    if (off=0) and (curr_blk=0) then begin
      {special case 2}
      curr_off := 1;
      prime_sieve_next := 2;
      exit;
    end;
    repeat
      while off < SIEVE_BLOCK do begin
        if not psieve^[off] then begin
          curr_off := off + 1;
          {If SIEVE_BLOCK2 is no power of 2 use * instead of shl!}
          {prime_sieve_next := curr_blk * SIEVE_BLOCK2 + curr_off + off;}
          prime_sieve_next := curr_blk shl 16 + curr_off + off;
          exit;
        end;
        inc(off);
      end;
      inc(curr_blk);
      off := 0;
      if curr_blk<SIEVE_MAXBLOCK then prime_sieve_nextblock(sieve)
      else begin
        prime_sieve_next := 1;
        exit;
      end;
    until false;
  end;
end;


begin
  {If assert fails adjust code in prime_sieve_next}
  assert(SIEVE_BLOCK2 = 1 shl 16, MPAF+'SIEVE_BLOCK2 = 1 shl 16');
end.

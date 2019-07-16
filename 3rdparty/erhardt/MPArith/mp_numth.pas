unit mp_numth;

{MPArith number theoretic functions}

interface

{$i STD.INC}

{$ifndef FPC}
  {$N+}
{$endif}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  MPArith number theoretic functions

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25S, FPC, VP, WDOSX

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] LibTomMath 0.30+ by Tom St Denis
                  [2] MPI by M.J. Fromberger
                  [3] Knuth, D.E.: The Art of computer programming. Vol 2
                      Seminumerical Algorithms, 3rd ed., 1998
                  [4] Forster, O.: Algorithmische Zahlentheorie, 1996
                  [5] (HAC) Menezes,A., von Oorschot,P., Vanstone, S: Handbook of
                      Applied Cryptography, 1996, www.cacr.math.uwaterloo.ca/hac
                  [6] R. P. Brent, Factor: an integer factorization program for
                      the IBM PC, Report TR-CS-89-23, October 1989, 7 pp.
                      http://maths-people.anu.edu.au/~brent/pub/pub117.html
                      http://maths-people.anu.edu.au/~brent/ftp/rpb117/rpb117.exe
                  [7] P. Ribenboim: The New Book of Prime Number Records, 3rd ed., 1995.
                  [8] Marcel Martin: NX - Numerics library of multiprecision
                      numbers for Delphi and Free Pascal, 2006-2009
                      www.ellipsa.eu/public/nx/index.html
                  [9] Wei Dai: Lucas Sequences in Cryptography,
                      http://www.weidai.com/lucas.html
                 [10] Crandall,R., C.Pomerance: Prime Numbers, A Computational
                      Perspective, 2nd ed., 2005
                 [11] Peter Luschny's Homepage of Factorial Algorithms,
                      http://www.luschny.de/math/factorial/FastFactorialFunctions.htm
                 [12] PARI/GP at http://pari.math.u-bordeaux.fr/
                 [15] The GNU Multiple Precision Arithmetic Library, http://gmplib.org/
                 [24] H. Cohen, A Course in Computational Algebraic Number Theory
                      4th printing, 2000
                 [25] IEEE P1363/Draft. Standard Specifications for Public Key Cryptography.
                      Annex A (informative). Number-Theoretic Background.
                 [26] A. Adler, J.E. Coury: The Theory of Numbers, 1995
                 [32] S.C. Lindhurst, Computing Roots in Finite Fields and Groups, with a
                      Jaunt through Sums of Digits, University of Wisconsin, Madison 1997
                      http://scott.lindhurst.com/papers/thesis.ps.gz
                 [40] H. Riesel, Prime Numbers and Computer Methods for Factorization,
                      Vol. 126 of Progress in Mathematics, Boston, 2nd ed. 1994.
                      Paperback reprint 2012 in Modern Birkh„user Classics Series.


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.0.01   09.08.04  W.Ehrhardt  Initial version: mp_addmod, mp_submod
 0.0.02   09.08.04  we          mp_mulmod, mp_sqrmod
 0.0.03   09.08.04  we          raw versions of mp_pollard_rho, mp_gcd
 0.0.04   09.08.04  we          mp_pollard_rho with Barrett reduction
 0.0.05   09.08.04  we          mp_gcd (binary version), mp_exteuclid (xgcd)
 0.0.06   10.08.04  we          mp_pollard_rho: Barrett after x*x+2
 0.0.07   11.08.04  we          mp_lcm
 0.0.08   11.08.04  we          mp_invmod using xgcd
 0.0.08   12.08.04  we          mp_invmodb using binary gcd
 0.0.09   12.08.04  we          fast_mp_invmod, inv_mod -> inv_mod_euclid, inv_modb -> inv_mod
 0.0.10   12.08.04  we          code cleanup
 0.0.11   13.08.04  we          s_mp_exptmod
 0.0.12   13.08.04  we          basic mp_exptmod
 0.0.13   14.08.04  we          WEXP_ constants
 0.0.14   14.08.04  we          mp_pollard_rho with cnt parameter
 0.0.15   17.08.04  we          mp_is_square
 0.0.16   17.08.04  we          mp_small_factor
 0.0.17   17.08.04  we          mp_small_factor with mod 30 increments
 0.0.18   18.08.04  we          mp_is_SPP
 0.0.19   22.08.04  we          small prime tables
 0.0.20   22.08.04  we          mp_miller_rabin
 0.0.21   23.08.04  we          mp_show_progress variable, mp_small_factor with f0
 0.0.22   23.08.04  we          Prime127 removed, reordering of functions
 0.0.23   24.08.04  we          s_mp_ppexpo, mp_pollard_pm1
 0.0.24   24.08.04  we          code cleanup, updated references

 0.3.00   26.08.04  we          "release" version 0.3
 0.3.01   26.08.04  we          mp_error, most functions now procedures
 0.3.02   26.08.04  we          mp_progress, mp_set_progress
 0.3.03   16.02.05  we          Montgomery Reduction
 0.3.04   16.02.05  we          temp. $R- for mp_montgomery_setup (for 8-bit)
 0.3.05   17.02.05  we          general mp_exptmod_win, removed s_mp_exptmod
 0.3.06   20.02.05  we          mp_jacobi
 0.3.07   20.02.05  we          mp_jacobi: fix up D6/FPC warnings
 0.3.08   27.03.05  we          Montgomery routines to mp_core
 0.3.09   29.03.05  we          mp_jacobi: non-recursive version
 0.3.10   08.05.05  we          mp_jacobi: clean up, remove old version
 0.3.11   16.05.05  we          FPC: $mode objfpc/$goto on
 0.3.12   10.08.05  we          mp_lucasv.., debug check in mp_xgcd
 0.3.13   15.08.05  we          bugfix mp_lucasv2p if v same var as p or q
 0.3.13   15.08.05  we          mp_williams_pp1: William's p+1 method
 0.3.15   16.08.05  we          mp_coshmult compatible with Forster
 0.3.16   17.08.05  we          s_mp_coshmult easy out for b=1, loop error check
 0.3.17   18.08.05  we          mp_williams_pp1: numtest parameter

 0.4.00   20.08.05  we          use mp_set_error
 0.4.01   20.08.05  we          longint bit count in mp_exptmod_win
 0.4.02   21.08.05  we          MPC_ArgCheck
 0.4.03   22.08.05  we          removed mp_argcheck, mp_errchk
 0.4.04   22.08.05  we          now functions: mp_is_SPP, mp_is_square
 0.4.05   22.08.05  we          mp_fermat
 0.4.06   22.08.05  we          pollard, williams, coshmult with init<i>
 0.4.07   23.08.05  we          mp_mersenne, MPC_HaltOnArgCheck
 0.4.08   23.08.05  we          new values for MaxFermat, MaxMersenne
 0.4.09   23.08.05  we          IsPrime15 (quick and dirty prime test)
 0.4.10   23.08.05  we          mp_exptmod_d
 0.4.11   24.08.05  we          Fibonacci numbers, mp_fib, mp_fib2
 0.4.12   25.08.05  we          Lucas numbers, mp_lucas, mp_lucas2
 0.4.13   25.08.05  we          mp_fib2 for n=0
 0.4.14   26.08.05  we          miller_rabin: generate a in the range 2<=a<n-1
 0.4.15   27.08.05  we          easy outs in mp_small_factor
 0.4.16   28.08.05  we          usage of exceptions implemented
 0.4.17   10.09.05  we          function IsMersennePrime
 0.4.18   12.09.05  we          mp_fact, MaxFact
 0.4.19   12.09.05  we          optimized 3-level mp_fact, FSIVec
 0.4.20   13.09.05  we          mp_fact optimization (balanced factors)
 0.4.22   16.09.05  we          changed mp_gcd easy out, comment in mp_xgcd
 0.4.23   17.09.05  we          mp_ecm_brent (Brent's ECM factoring method)
 0.4.24   18.09.05  we          mp_ecpwr (faster version of mp_ecpwr2)
 0.4.25   20.09.05  we          mp_lucasv2p without goto, changed TProgressProc
 0.4.25   22.09.05  we          $argcheck for pointers (mp_lucasv2p, mp_xgcd)
 0.4.26   25.09.05  we          Bugfix: k now longint in IsMersennePrime
 0.4.27   26.09.05  we          $argcheck for identical args in mp_fib2,mp_lucas2

 0.5.00   29.09.05  we          'internal' functions moved to end of interface
 0.5.01   30.09.05  we          IsPrime16, pbits16, pmask16
 0.5.02   01.10.05  we          easy out for IsMersennePrime
 0.5.03   02.10.05  we          changed 'SPP' to more conventional 'spsp'
 0.5.04   03.10.05  we          function nextprime32
 0.5.05   03.10.05  we          function prevprime32
 0.5.06   09.10.05  we          next/prev prime residues classes module via $ifdef
 0.5.07   09.10.05  we          mp_is_spsp_d, mp_is_pprime, mp_small_factor with fmax
 0.5.08   09.10.05  we          small changes in mp_miller_rabin
 0.5.08   13.10.05  we          mp_is_slpsp
 0.5.09   16.10.05  we          mp_is_pprime with BPSW probable prime
 0.5.10   14.11.05  we          mp_nextprime/mp_prevprime
 0.5.11   15.11.05  we          mp_rand_prime
 0.5.12   18.11.05  we          mp_rand_prime improved for safe primes
 0.5.13   20.11.05  we          TPrimeType = (pt_normal, pt_3mod4, pt_safe)
 0.5.14   20.11.05  we          mp_is_pprime_ex, mp_rand_prime with mp_is_pprime_ex
 0.5.15   21.11.05  we          mp_rand_prime with mp_rand_bits
 0.5.16   22.11.05  we          New name: mp_isMersennePrime, mp_prng removed from uses

 0.6.00   30.12.05  we          mp_count_bits/CountBits32 renamed to mp_bitsize/bitsize32
 0.6.01   30.12.05  we          MP_8BIT removed
 0.6.02   05.01.06  we          MaxMersenne = DIGIT_BIT*MAXDigits-1;
 0.6.03   10.01.06  we          mp_fact using Recursive Split

 0.7.00   28.04.06  we          mp_xgcd: u2,v2,t2 suppressed
 0.7.01   04.08.06  we          TRedType MR_Reduce2k, mp_reduce_2k in mp_exptmod_win
 0.7.03   05.08.06  we          mp_show_progress initialized with false
 0.7.04   08.08.06  we          mp_makeodd in mp_miller_rabin, mp_is_slpsp, mp_jacobi
 0.7.05   08.08.06  we          mp_sqrtmod initial version
 0.7.06   09.08.06  we          mp_sqrtmod arg/init checks
 0.7.07   09.08.06  we          mp_lucasvmod
 0.7.08   09.08.06  we          mp_sqrtmod18_lucas, TPrimeType with pt_1mod8
 0.7.09   10.08.06  we          mp_sqrtmodp2, mp_sqrtmethod
 0.7.10   11.08.06  we          mp_miller_rabin and mp_is_spsp modified
 0.7.11   11.08.06  we          Avoid FPC warnings: nextprime32/prevprime32
 0.7.12   11.08.06  we          experimental mp_sqrtmod18_fp
 0.7.13   11.08.06  we          check prime/q-residue after 100 trials in mp_sqrtmod18_shanks/fp
 0.7.14   15.08.06  we          small tweak in mp_isMersennePrime
 0.7.15   15.08.06  we          $ifdef MPC_sqrtmod_fp
 0.7.16   16.08.06  we          bigprime_pm1: big prime stage for Pollard p-1
 0.7.17   17.08.06  we          bigprime_pm1: Barret and bitsize dependent constants
 0.7.17   20.08.06  we          const mp_max_small / mp_small_factor
 0.7.18   21.08.06  we          mp_ecm_factor (first complete version using 7KB stack)
 0.7.19   22.08.06  we          mp_ecm_factor (dynamic vectors, progress, intermeditate GCDs)
 0.7.20   23.08.06  we          mp_ecm_factor,mp_ecm_brent: check n>1, @n<>@f
 0.7.21   23.08.06  we          s_mp_ppexpo with mp_mul_w
 0.7.22   23.08.06  we          mp_ecm_factor with Barrett: speed more than doubled
 0.7.23   24.08.06  we          mp_ecm_factor: check seed>5 in ecm_setup, more comments
 0.7.24   26.08.06  we          nextprime32_array
 0.7.25   27.08.06  we          Prime table in mp_ecm_factor $ifdef MPC_ECM_Primetable
 0.7.26   28.08.06  we          basic mp_is_power
 0.7.27   07.09.06  we          mp_is_square2
 0.7.28   07.09.06  we          nextprime32_array with parameter cmax
 0.7.29   07.09.06  we          mp_ecm_factor with parameters C1,phase; mp_ecm_simple
 0.7.30   07.09.06  we          function ProgressAssigned

 0.8.00   19.09.06  we          mp_primorial, MaxPrimorial
 0.8.01   19.09.06  we          mp_isMersennePrime: search small factor before Lucas/Lehmer
 0.8.02   23.09.06  we          MaxFibonacci, MaxLucas
 0.8.03   25.09.06  we          mp_kronecker, mp_jacobi uses kron_intern

 0.9.00   25.12.06  we          minor changes in mp_isMersennePrime
 0.9.01   28.12.06  we          Bugfix mp_is_power for a=-4
 0.9.02   30.12.06  we          FindFirst/NextPrime32, rewrite nextprime32_array
 0.9.03   01.01.07  we          bigalloc renamed to IAlloc
 0.9.04   01.01.07  we          mp_primorial with mp_prod_int
 0.9.05   01.03.07  we          mp_primorial with mp_primor_cutoff

 1.0.00   07.05.07  we          mp_exptmod returns 0 for c=1, error for c<=0
 1.0.01   08.05.07  we          Removed Fp[sqrt(D)] arithmetic for MPC_sqrtmod_fp
 1.0.02   09.05.07  we          Removed mp_ecpwr2
 1.0.03   13.05.07  we          Corrected some exception strings

 1.1.00   19.05.07  we          break mp_gcd loop if u=v
 1.1.01   22.05.07  we          Binary Kronecker symbol algorithm
 1.1.02   22.05.07  we          Binary Kronecker: skip operations after mp_sub
 1.1.03   23.05.07  we          Binary Kronecker: s_mp_sub, break if x=y
 1.1.04   31.05.07  we          mp_is_power: check mod 4 before mp_expt_int
 1.1.05   27.06.07  we          mp_gcd_initial_mod in mp_gcd
 1.1.06   01.07.07  we          removed useless first iteration in mp_fact
                                (thanks to Marcel Martin)
 1.1.07   02.07.07  we          mp_miller_rabin uses IsPrime16 if possible
 1.1.08   09.07.07  we          mp_gcd_initial_mod renamed to mp_initial_mod
                                and also used in kron_intern
 1.1.09   12.07.07  we          mp_pell and mp_pell4
 1.1.10   14.07.07  we          s_mp_lucasvmod1 (used by mp_is_slpsp)
 1.1.11   14.07.07  we          mp_sqrtmod14_mueller
 1.1.12   14.07.07  we          Mueller sqrtmod: Bugfix use 1/t mod p
 1.1.13   15.07.07  we          Mueller sqrtmod: case t=1 handled in Jacobi loop
 1.1.14   15.07.07  we          Mueller sqrtmod: don't compute PP in Jacobi loop
 1.1.15   21.07.07  we          mp_xgcd: easy out for trivial cases
 1.1.17   30.07.07  we          mp_sqrtmodpk
 1.1.17   31.07.07  we          mp_sqrtmod14_mueller for mp_sqrtmethod=3 and
                                in auto mode instead of Lucas
 1.1.18   01.08.07  we          fast_mp_invmod: no init checks, y removed
 1.1.19   04.08.07  we          mp_sqrtmod2k
 1.1.20   05.08.07  we          special cases k=1,2 of mp_sqrtmod2k
                                easy out in mp_exptmod if b=1
 1.1.21   05.08.07  we          mp_sqrtmod_ex with optional Jacobi check
 1.1.22   05.08.07  we          mp_sqrtmodpq, removed mp_invmod_euclid

 1.2.00   16.08.07  we          mp_lcm: special cases, divide larger arg by gcd
 1.2.01   17.08.07  we          mp_is_primepower
 1.2.02   21.08.07  we          jacobi32
 1.2.03   23.08.07  we          mp_jacobi_lm
 1.2.04   24.08.07  we          mp_jacobi_lm used in: mp_is_slpsp, mp_sqrtmod18_shanks_intern
 1.2.05   26.08.07  we          improved error handling in mp_sqrtmod18_* routines
 1.2.06   30.08.07  we          t=1 in mp_sqrtmod14_mueller, improved jacobi32
 1.2.07   02.09.07  we          updated mp_primor_cutoff values
 1.2.08   03.09.07  we          mp_OddProd, mp_dfact, mp_fact with local mp stack
 1.2.09   04.09.07  we          small tables for mp_lucas and mp_fib
 1.2.10   04.09.07  we          Bugfix mp_fib2 for n=k*$10000000, k=1..7
 1.2.11   06.09.07  we          mp_binomial, s_mp_binom_w
 1.2.12   07.09.07  we          mp_binomial for negative arguments
 1.2.13   08.09.07  we          bugfix mp_binomial for small k and n>65535
 1.2.14   08.09.07  we          s_mp_binom_l, case k=1 in mp_binomial
 1.2.15   09.09.07  we          mp_binomial: init X[s] when needed
 1.2.16   10.09.07  we          mp_sqrtmod916_kong
 1.2.17   10.09.07  we          avoid warning in mp_miller_rabin if MP_16BIT
 1.2.18   11.09.07  we          mp_jacobi_ml
 1.2.19   12.09.07  we          Bugfix mp_is_primepower for a=2^1
 1.2.20   17.09.07  we          mp_is_power uses residues mod 8

 1.3.00   26.10.07  we          s_mp_sqrtmod2k for odd a, mp_sqrtmod2k for a>=0
 1.3.01   12.11.07  we          Fix memory leak(s) if MPC_HaltOnError is not defined
 1.3.02   10.12.07  we          improved mp_is_square2

 1.5.00   13.01.08  we          mp_sqrtmod2k: bugfix due to error in Adler/Coury [26]

 1.6.00   24.05.08  we          mp_crt_solve, mp_crt_setup, mp_crt_single
 1.6.01   24.05.08  we          functions mp_invmodf, mp_crt_setupf
 1.6.02   04.06.08  we          mp_is_square2: psqrt^ := a if a=0 or 1
 1.6.03   04.06.08  we          mp_cornacchia
 1.6.04   05.06.08  we          mp_cornacchia4
 1.6.05   06.06.08  we          function mp_gcd1
 1.6.06   08.06.08  we          Cosmetic changes / corrected exception strings
 1.6.07   11.06.08  we          Arg check in mp_binomial

 1.7.00   27.07.08  we          Improved trial factors in mp_isMersennePrime
 1.7.01   17.09.08  we          use mp_is_longint after mp_initial_mod in mp_gcd
 1.7.02   17.09.08  we          mp_gcd_euclid with mp_mod, more uses of mp_is_longint
 1.7.03   20.09.08  we          mp_xgcd with Lehmer, mp_gcd_ml (modified Lehmer)
 1.7.04   21.09.08  we          $ifdef MPC_UseGCD32 in mp_gcd
 1.7.05   22.09.08  we          mp_xgcd_bin
 1.7.06   22.09.08  we          mp_invmodf via xgcd, removed fast_mp_invmod
 1.7.07   22.09.08  we          extended range for Lehmer routine
 1.7.08   23.09.08  we          improved inner loop of (binary) mp_gcd
 1.7.09   24.09.08  we          use mp_sign for trivial cases in mp_xgcd_bin
 1.7.10   24.09.08  we          mp_is_pprime_ex skips 1-digit test, fixed mp_is_square2
 1.7.11   24.09.08  we          mp_invmodf with invmod32
 1.7.12   02.10.08  we          bugfix in mp_gcd $ifdef MPC_UseGCD32

 1.8.00   09.10.08  we          use mp_set1; mp_is_longint in mp_miller_rabin
 1.8.01   21.10.08  we          exptmod32, s_mp_is_pth_power, improved mp_is_power
 1.8.02   22.10.08  we          changed mp_is_primepower
 1.8.03   25.10.08  we          mp_is_power with s_mp_n_root2
 1.8.04   26.10.08  we          Check a mod k^3 and a mod q^2 in s_mp_is_pth_power
 1.8.05   27.10.08  we          improved mp_is_power
 1.8.06   01.11.08  we          mp_xlcm, mp_rnr
 1.8.07   02.11.08  we          sieve test loop in s_mp_is_pth_power, new: mp_is_pth_power

 1.9.00   06.12.08  we          mp_rnr2
 1.9.01   23.12.08  we          s_mp_mca_alg1816
 1.9.02   25.12.08  we          improved mp_small_factor
 1.9.03   29.12.08  we          uses mp_prime
 1.9.04   03.01.09  we          mp_is_square2, s_mp_pell4 with s_mp_sqrtrem
 1.9.05   04.01.09  we          mp_is_square2: s_mp_sqrtrem only (removed trick branch)
                                removed mp_ecm_brent, mp_ecpwr
 1.9.06   04.01.09  we          32 bit first/next prime routines to mp_prime
 1.9.07   06.01.09  we          minor changes in mp_isMersennePrime
 1.9.08   06.01.09  we          explicit range check in mp_fermat

 1.10.00  21.01.09  we          changes related to (s)mp_divrem
 1.10.01  22.01.09  we          mp_cbrtmod
 1.10.02  25.01.09  we          s_mp_is_pprime_ex
 1.10.03  28.01.09  we          improved mp_is_square2
 1.10.04  02.02.09  we          mp_sqrtmodpk: check invmod, store result only if Err=0
 1.10.05  03.02.09  we          s_mp_sqrtmodpk with red parameter, mp_sqrtmodpk with red=true
 1.10.06  03.02.09  we          mp_cbrtmodpk, s_mp_cbrtmodpk, mp_cbrtmod3k
 1.10.07  04.02.09  we          sqrt/cbrtmodpk functions: handle a mod p = 0
 1.10.08  04.02.09  we          mp_sqrtmodpk handles p=2
 1.10.09  07.02.09  we          mp_cbrtmodpq
 1.10.10  07.02.09  we          s_mp_is_cubres
 1.10.11  08.02.09  we          mp_cbrtmod_ex
 1.10.12  12.02.09  we          improved mp_cbrtmod_ex
 1.10.13  19.02.09  we          improved s_mp_is_cubres
 1.10.14  27.02.09  we          mp_isMersennePrime: check for spurious factor

 1.12.00  20.06.09  we          mp_is_power_max
 1.12.01  21.06.09  we          mp_provable_prime
 1.12.02  05.07.09  we          new: s_mp_npfirst, s_mp_npnext, updated: mp_nextprime, s_mp_is_pprime_ex
 1.12.04  29.07.09  we          Increased trace level in mp_provable_prime

 1.13.00  15.08.09  we          mp_crt_single/solve: improved and check n>0
 1.13.01  23.08.09  we          allow p=q in mp_sqrtmodpq/mp_cbrtmodpq

 1.16.00  04.06.10  we          mp_4sq, mp_4sq_sa, mp_4sq_sd

 1.17.00  30.12.10  we          mp_binomial for k < 0

 1.19.00  03.11.11  we          Adjust MaxFermat, MaxFact for MAXDigits=32000 with MP_32BIT

 1.20.00  14.01.12  we          BIT64: exptmod32
 1.20.01  21.01.12  we          Adjusted MaxPrimorial, MaxPrimorial, NMAX in mp_primorial

 1.21.00  15.07.12  we          mp_isconsole: boolean [used for write('.')]
 1.21.01  18.07.12  we          mp_catalan, MaxCatalan
 1.21.02  19.07.12  we          improved mp_small_factor
 1.21.03  20.07.12  we          mp_poch
 1.21.04  20.07.12  we          mp_perm
 1.21.05  21.07.12  we          mp_poch and mp_perm use mp_product
 1.21.06  25.07.12  we          exptmod32/jacobi32 moved to mp_base
 1.21.07  27.07.12  we          mp_is_spsp returns false for n<=1
 1.21.08  29.07.12  we          mp_is_square2 uses is_square32ex
 1.21.09  29.07.12  we          mp_squfof
 1.21.10  30.07.12  we          mp_sigmak
 1.21.11  31.07.12  we          mp_sigmak with geometric sum to avoid overflow

 1.22.00  10.08.12  we          mp_nextprime_ex
 1.22.01  11.08.12  we          mp_small_factor with inline IsPrime16
 1.22.02  12.08.12  we          mp_pollard_pm1 : max. bound = 65000, fix while condition bug
                                mp_williams_pp1: max. bound = 65000
 1.22.03  14.08.12  we          Borwein/Schoenhage s_mp_borsch_fact
 1.22.04  15.08.12  we          mp_pollard_brent/ex
 1.22.05  16.08.12  we          Changed mp_lucasv2p; new: mp_lucasuv, mp_lucasu
 1.22.06  18.08.12  we          Removed mp_coshmult
 1.22.07  19.08.12  we          mp_squad_mod
 1.22.08  21.08.12  we          mp_cornacchia_ex
 1.22.09  27.08.12  we          mp_qnr
 1.22.10  01.09.12  we          special treatment for d=12 in s_mp_pell4
 1.22.11  06.09.12  we          mp_rqffu and s_mp_rqffu
 1.22.12  08.09.12  we          mp_qnr for all n > 0
 1.22.13  11.09.12  we          mp_powerd
 1.22.14  12.09.12  we          mp_ecm_factor uses Primes16

 1.23.00  16.09.12  we          fix s1 in mp_small_factor
 1.23.01  24.09.12  we          move s_mp_mca_alg1816 to mp_rsa
 1.23.02  24.09.12  we          s_mp_nextprime_sieve
 1.23.03  24.09.12  we          basic modular arithmetic, GCD/LCM moved to mp_modul
 1.23.04  24.09.12  we          factorization routines moved to mp_pfu
 1.23.05  28.09.12  we          mp_cornacchia_ex renamed to mp_cornacchia_pk
 1.23.06  29.09.12  we          s_mp_cornacchia_ex, mp_cornacchia_pq
 1.23.07  30.09.12  we          improved mp_is_power
 1.23.08  02.10.12  we          BIT: 16fix case in s_mp_is_pth_power with p > 2^16
 1.23.09  03.10.12  we          mp_val, mp_valrem, s_mp_valrem_d
 1.23.10  03.10.12  we          mp_is_power with mp_valrem
 1.23.11  04.10.12  we          Fix bug for test d=5,12 in s_mp_pell4
 1.23.12  14.10.12  we          mp_is_power with flexible trial division

 1.24.00  17.12.12  we          Some word types changed to longint
 1.24.01  19.12.12  we          mp_fact with longint arguments, s_mp_borsch_fact uses sieve
 1.24.02  20.12.12  we          MaxFact, MaxFermat computed in initialization
 1.24.03  20.12.12  we          mp_binomial with prime_sieve for 32+ bits
 1.24.04  20.12.12  we          mp_dfact and mp_OddProd with longint arguments
 1.24.05  21.12.12  we          s_mp_gbinom (used in mp_binomial and mp_catalan)
 1.24.06  03.01.13  we          Improved mp_primorial
 1.24.07  04.01.13  we          Improved s_mp_sqrtmod2k; mp_sqrtmod2k: No error if a is a multiple of 2^k
 1.24.08  04.01.13  we          Handle k=0 in s_mp_cbrtmodpk, error if k<0

 1.27.00  04.09.13  we          Allow a < 0 in mp_sqrtmod2k
 1.27.01  08.02.14  we          mp_val/rem: fix for a=0, changed iterative logic

 1.29.00  15.07.14  we          s_mp_nextprime_sieve can return safe primes
 1.29.01  16.07.14  we          mp_safeprime
 1.29.02  16.07.14  we          MPC_SieveNextPrime: use sieve for mp_next_prime
 1.29.03  17.07.14  we          s_mp_is_psp2
 1.29.04  17.07.14  we          s_mp_nextprime_sieve with s_mp_is_psp2
 1.29.05  18.07.14  we          mp_is_spsp with Barret reduction if s >= SBMin
 1.29.06  18.07.14  we          mp_rand_prime with mp_safeprime for pt_safe type

 1.30.00  29.07.14  we          Improved s_mp_sqrtmodpk for k a power of 2
 1.30.01  11.09.14  we          Improved s_mp_cbrtmodpk for k a power of 2
 1.30.02  20.09.14  we          mp_is_power_modpk
 1.30.03  21.09.14  we          s_mp_nroot_modp
 1.30.04  21.09.14  we          s_mp_nroot_modpk
 1.30.05  22.09.14  we          s_mp_nroot_modpq
 1.30.06  25.09.14  we          Fixed mp_is_power_modpk with n: mp_int
 1.30.07  28.09.14  we          mp_safeprime_ex
 1.30.08  29.09.14  we          Improved mp_perm
 1.30.09  29.09.14  we          Check/handle n<0 in mp_fact

 1.31.00  17.10.14  we          small change in s_mp_is_pth_power
 1.31.01  05.11.14  we          Fix/extend mp_dfact for n<0
 1.31.02  22.11.14  we          Extended s_mp_nroot_modp

 1.32.00  03.03.15  we          Exception strings corrected in mp_cbrtmod3k, mp_cbrtmod_ex

 1.33.00  14.09.16  we          Allow n > p-1 in s_mp_nroot_modp

 1.34.00  09.12.16  we          mp_is_pcarmichael/ex

 1.35.00  25.07.17  we          mp_fact: check n > MaxFact
 1.35.01  26.07.17  we          mp_fact_swing
 1.35.02  26.07.17  we          mp_swing

 1.36.00  10.09.17  we          mp_is_psp
 1.36.01  11.09.17  we          mp_is_psp_d
 1.36.02  12.09.17  we          mp_nextpower/_ex
 1.36.03  18.09.17  we          mp_digitsum
 1.36.04  21.09.17  we          s_mp_is_lpsp with Selfridge A parameters
 1.36.05  22.09.17  we          Rename existing strong Lucas test to mp_is_slpsp_alt
 1.36.06  22.09.17  we          mp_is_lpsp, mp_is_slpsp with s_mp_is_lpsp

 1.37.00  15.10.17  we          mp_reverse

**************************************************************************)

(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - LibTomMath 0.30+ by Tom St Denis
   - MPI 1.8.6 by Michael J. Fromberger
   - NX V0.18 and V0.9+ by Marcel Martin
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2017 Wolfgang Ehrhardt

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


const
  mp_max_small = mp_digit(MP_DIGIT_MAX and $7FFF); {max. tested small factor}

var
  MaxFact: longint;           {MaxFact depends on DIGIT_BIT and MAXDigits. It  }
                              {is an approx. value for max. valid mp_fact args.}
                              {Since no checks are performed some larger values}
                              {might work. With 16 bit real mode even smaller  }
                              {values may result in "Out of memory" errors.    }
var
  MaxFermat: integer;         {Maximum Fermat number with this implementation  }

var
  mp_sqrtmethod   : integer;  {Selects algorithm for sqrt mod p, p=1 mod 8.}
                              {0:Auto, 1:Shanks, 2:Lucas, 3:Mueller.}
                              {Default=0: auto select between Shanks/Mueller.}

var
  mp_primor_cutoff: longint;  {Cutoff point for recursive product in mp_primorial,}
                              {initialization will set rough value. For better}
                              {results use t_tunepr}

type
  TPrimeType    = (pt_normal, pt_3mod4, pt_safe, pt_1mod8);
                              {pt_normal: prime without special features.}
                              { pt_3mod4: prime p with p modulo 4 = 3.   }
                              { pt_1mod8: prime p with p modulo 8 = 1.   }
                              {  pt_safe: safe prime, ie p=2q+1, q prime.}


const
  MaxMersenne  = DIGIT_BIT*MAXDigits-1;      {Approx. maximum Mersenne number}
  MaxPrimorial = trunc(MaxMersenne*0.6932);  {Approx. maximum arg for mp_primorial, Note: for real}
                                             {mode "out of memory" will occur for smaller values!}
  MaxFibonacci = trunc(MaxMersenne*1.4404);  {Approx. maximum arg for mp_fib}
  MaxLucas     = MaxFibonacci-2;             {Approx. maximum arg for mp_lucas}
  MaxCatalan   = trunc(MaxMersenne*0.4999);  {Approx. maximum arg for mp_catalan}


procedure mp_4sq(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}

procedure mp_4sq_sa(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2, a<=b<=c<=d.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}

procedure mp_4sq_sd(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2, a>=b>=c>=d.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}

procedure mp_binomial(n,k: longint; var a: mp_int);
  {-Calculate the binomial coefficient a = (n choose k)}

procedure mp_catalan(n: longint; var c: mp_int);
  {-Return the n'th Catalan number C_n = binomial(2n,n)/(n+1), 0 if n < 0}

procedure mp_cbrtmod(const a, p: mp_int; var x: mp_int; var Err: integer);
  {-Compute a cube root x of a with x^3 = a mod p, p prime.  }
  { Err in [1..3], if a is not a 3rd power or p is not prime.}

procedure mp_cbrtmod3k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod 3^k. Error codes: 1<Err<4 code from}
  { mp_cbrtmod, Err=4 if p=3, Err=5: no inverse mod 3^(2^i), Err=6: (internal d<>0)}

procedure mp_cbrtmodpk(const a,p: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod p^k, p prime. b will be  }
  { reduced mod p^k. Error codes: 1<Err<4: from mp_cbrtmod, Err=4 if p=3,}
  { Err=5: no inverse mod p^(2^i), Err=6: internal check for p=3 failed. }

procedure mp_cbrtmodpq(const a,p,q: mp_int; var x: mp_int; var Err: integer);
  {-Compute a cube root x of a mod (pq); p,q primes. If p=q, mp_cbrtmodpk is}
  { used. For p<>q: Err=8 if gcd(p,q)<>1, otherwise error code from mp_cbrtmod.}

procedure mp_cbrtmod_ex(const a, p: mp_int; var x: mp_int; pru: pmp_int; var Err: integer);
  {-Compute a cube root x of a with x^3 = a mod p, p prime. If pru<>nil, pru^ }
  { is set to a 3rd root of unity (or 1 if x is unique), the other roots are  }
  { x*pru^, x*pru^*pru^. Err <> 0, if a is not a 3rd power or p is not prime. }

procedure mp_cornacchia(const d,p: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve x^2 + d*y^2 = p, p prime, 0<d<p with Cornacchia's algorithm. Check}
  { @x<>@y, but not p prime. Err<>0 if no solution exist. x >= y if d=1.}

procedure mp_cornacchia_pk(const a,b,p: mp_int; k: longint; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = p^k, p prime; k,a,b > 0, a + b < p^k, (a,p)=1 with Cornacchia's  }
  { algorithm. Check @x<>@y, but not p prime. Err <> 0 if no solution exist. x >= y if a=b.}

procedure mp_cornacchia_pq(const a,b,p,q: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = p*q, p,q prime; a,b > 0, a + b < p*q, (a,p*q)=1 with Cornacchia's  }
  { algorithm. Check @x<>@y, but not p,q prime. Err <> 0 if no solution exist. x >= y if a=b.}

procedure mp_cornacchia4(const d,p: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve x^2 + |d|*y^2 = 4p, p prime, -4p<d<0 with the modified Cornacchia}
  { algorithm. Check @x<>@y, but not p prime. Err<>0 if no solution exist.}

procedure mp_crt_setup(n: integer; const m: array of mp_int; var c: array of mp_int);
  {-Calculate CRT coefficients c[i] for pairwise co-prime moduli m[i], i=0..n-1.}

function  mp_crt_setupf(n: integer; const m: array of mp_int; var c: array of mp_int): boolean;
  {-Calculate CRT coefficients c[i] for pairwise co-prime moduli m[i], i=0..n-1.}
  { Return true if c[i] are successfully calculated.}

procedure mp_crt_single(n: integer; const m,v: array of mp_int; var x: mp_int);
  {-Calculate x with x mod m[i] = v[i], pairwise co-prime moduli m[i], i=0..n-1.}
  { Use mp_crt_setup/calc if more than one system shall be solved}

procedure mp_crt_solve(n: integer; const m,c,v: array of mp_int; var x: mp_int);
  {-Calculate x with x mod m[i] = v[i], pairwise co-prime moduli m[i], i=0..n-1.}
  { Coefficients c[i] must be precalculated with mp_crt_setup}

function  mp_digitsum(const a: mp_int; radix: word): longint;
  {-Return the digitsum of a in base radix (2..MAXRadix)}

procedure mp_dfact(n: longint; var a: mp_int);
  {-Calculate double factorial a = n!!, a=n*(n-2)..., lowest term 1 or 2; a=0 for n < -3}

procedure mp_fact(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) using Recursive Split or Borwein/Schoenhage}
  { prime factorization method, error if n > MaxFact}

procedure mp_fact_swing(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) using Luschny's prime swing method}

procedure mp_fermat(n: word; var fn: mp_int);
  {-Return nth Fermat number, fn = 2^(2^n)+1 (MP_RANGE error for n>MaxFermat)}

procedure mp_fib(n: longint; var fn: mp_int);
  {-Calculate Fibonacci number fn=fib(n), fib(-n)=(-1)^(n-1)*fib(n)}

procedure mp_fib2(n: longint; var fn,f1: mp_int);
  {-Calculate two Fibonacci numbers fn=fib(n), f1=fib(n-1), n>=0}

function  mp_isMersennePrime(p: longint): boolean;
  {-Lucas-Lehmer test for Mersenne number m=2^p-1, HAC 4.37}

function  mp_is_lpsp(const a: mp_int): boolean;
  {-Lucas probable prime test for n using Selfridge A method}

function  mp_is_pcarmichael(const a: mp_int): boolean;
  {-Test if a is is a Carmichael number with high probability, 100 tests, prime check}

function  mp_is_pcarmichael_ex(const a: mp_int; numtests: integer; chkprime: boolean): boolean;
  {-Test if a is a Carmichael number with high probability, perform numtests}
  { tests (100 if numtests < 1) and optionally checks if a is prime.}

procedure mp_is_power(const a: mp_int; var b: mp_int; var p: longint);
  {-Calculate smallest prime p with a=b^p; p=1,b=a if a is no power}

procedure mp_is_power_max(const a: mp_int; var b: mp_int; var k: longint);
  {-Calculate largest k with a=b^k; k=1,b=a if a is no power}

function  mp_is_power_modpk(a,n,p: mp_int; k: longint): boolean;
  {-Return true if a is a nth power residue mod p^k;}
  { p prime (not checked); gcd(a,p)=1; n, k > 0.    }

function  mp_is_pprime(const a: mp_int): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32)}

function  mp_is_pprime_ex(const a: mp_int; smax: mp_digit): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32); trial division up to smax}

function  mp_is_primepower(const a: mp_int; var b: mp_int; var k: longint): boolean;
  {-Return true if a=b^k, b prime, k>1, otherwise false and a=b^k, k=1 if no power}

function  mp_is_pth_power(const a: mp_int; p: longint; var r: mp_int): boolean;
  {-Return true if a is pth power, a>0, p prime. If true, calculate r with a=r^p}

function  mp_is_psp(const n,a: mp_int): boolean;
  {-Probable prime test of n to base a > 1, true if a^(n-1) mod n = 1}

function  mp_is_psp_d(const n: mp_int; a: mp_digit): boolean;
  {-Probable prime test of n to base a > 1, true if a^(n-1) mod n = 1}

function  mp_is_slpsp_alt(const a: mp_int): boolean;
  {-Strong Lucas probable prime test for a. Lucas test is }
  { done for the first p=2k+1 with mp_jacobi(p^2-4,a) = -1}

function  mp_is_slpsp(const a: mp_int): boolean;
  {-Lucas (strong) probable prime test for n using Selfridge A method}

function  mp_is_spsp(const n,a: mp_int): boolean;
  {-Strong probable prime test of n to base a > 1 from HAC p. 139 Alg.4.24}

function  mp_is_spsp_d(const n: mp_int; a: mp_digit): boolean;
  {-Strong probable prime test of n to mp_digit base a > 1 from HAC p. 139 Alg.4.24}

function  mp_is_square(const a: mp_int): boolean;
  {-Test if a is square}

function  mp_is_square2(const a: mp_int; psqrt: pmp_int): boolean;
  {-Test if a is square, return sqrt(a) if a is a square and psqrt<>nil}

procedure mp_lucas(k: longint; var lk: mp_int);
  {-Calculate Lucas number lk=luc(k), luc(-k)=(-1)^k*luc(k)}

procedure mp_lucas2(k: longint; var lk,l1: mp_int);
  {-Calculate two Lucas numbers lk=luc(k), l1=luc(k-1), k>=0}

procedure mp_lucasv(const p,q: mp_int; k: longint; var v: mp_int);
  {-Calculate v[k] of Lucas V sequence for p,q, p^2-4q <>0, k>=0}

procedure mp_lucasuv(const p,q: mp_int; k: longint; var uk,vk: mp_int);
  {-Calculate u[k], v[k] of Lucas sequence for p,q; p^2-4q<>0, k>=0}

procedure mp_lucasu(const p,q: mp_int; k: longint; var u: mp_int);
  {-Calculate u[k] of Lucas sequence for p,q, p^2-4q <>0, k>=0}

procedure mp_lucasvmod(const p,q,n,k: mp_int; var vk: mp_int);
  {-Calculate v[k] mod n of Lucas V sequence for p,q.}
  { Ranges n>1, 0 <= p,q,k < n (checked if MPC_ArgCheck).}
  { Note: p^2-4q<>0 is not tested, no proper Lucas sequence!!}

procedure mp_lucasv2(const p,q: mp_int; k: longint; var v,w: mp_int);
  {-Calculate v=v[k],w=v[k+1] of Lucas V sequence for p,q, p^2-4q<>0, k>=0}

procedure mp_lucasv2p(const p,q: mp_int; k: longint; var vk: mp_int; p2: pmp_int; uk: boolean);
  {-Calculate v[k] of Lucas V sequence for p,q; p^2-4q<>0, k>=0;}
  { if p2<>nil, p2^ will be set to v[k+1] or u[k] if uk is true.}

procedure mp_mersenne(n: longint; var mn: mp_int);
  {-Return nth Mersenne number, mn = 2^n-1, MP_RANGE err for n>MaxMersenne}

procedure mp_miller_rabin(const n: mp_int; t: word; var prime: boolean);
  {-Miller-Rabin test of n, security parameter t, from HAC p. 139 Alg.4.24}
  { if t<=0, calculate t from bit_count with prob of failure < 2^-96}

procedure mp_nextpower(var a: mp_int);
  {-Return smallest perfect power >= a}

procedure mp_nextpower_ex(const a: mp_int; var b,c: mp_int; var p: longint);
  {-Return smallest perfect power b := c^p >= a}

procedure mp_nextprime(var n: mp_int);
  {-Next prime >= n, 2 if n<=2}

procedure mp_nextprime_ex(var n: mp_int; smax: mp_digit);
  {-Next prime >= n, 2 if n<=2; trial division from 3 to smax}

procedure mp_OddProd(a,b: longint; var p: mp_int);
  {-Calculate p=prod(2*i+1),i=a+1...b;  p=1 if a>=b}

procedure mp_pell(const d: mp_int; var x,y: mp_int);
  {-Calculate the smallest non-trivial solution of the Pell equation}
  { x^2 - d*y^2 = 1; error if d<2, if d is a square, or if @x = @y.}

procedure mp_pell4(const d: mp_int; var x,y: mp_int; var r: integer);
  {-Calculate the smallest non-trivial solution of  x^2 - d*y^2 = r,}
  { r in [-4,+4,-1,+1]. Uses continued fraction expansion of sqrt(d)}
  { from Forster [4]. Error if d<2, if d is a square, or if @x = @y.}

procedure mp_perm(n,r: longint; var a: mp_int);
  {-Compute a=n!/(n-r)!, the number of permutations of n distinct objects}
  { taken r at a time, n,r >= 0.}

procedure mp_poch(n,k: longint; var a: mp_int);
  {-Return the Pochhammer symbol a = (n)_k = n*(n+1)*...(n+k-1), a=1 if k<1}
  { (n)_k is sometimes called "rising factorial" or "Pochhammer function". }

procedure mp_powerd(const p,q,d: mp_int; n: longint; var x,y: mp_int);
  {-Calculate x + y*sqrt(d) = (p + q*sqrt(d))^n, n >= 0}

procedure mp_prevprime(var n: mp_int);
  {-Previous prime <= n, 0 if n<2}

procedure mp_primorial(n: longint; var a: mp_int);
  {-Primorial of n;  a = n# = product of primes <= n.}

procedure mp_provable_prime(bits: longint; var p: mp_int);
  {-Generate a random provable prime p with bitsize bits using Maurer's algorithm}

function  mp_qnr(const n: mp_int): longint;
  {-Return a small quadratic nonresidue for n > 0, -1 if error or no QNR is}
  { found: the smallest prime with kronecker(p|n) = -1 is returned. If n is}
  { a prime, the result is the least positive QNR, for composite n this may}
  { differ from the least QNR: e.g. mp_qnr(15) returns 7 and not 2.}

procedure mp_rand_prime(bitsize: longint; pt: TPrimeType; var p: mp_int);
  {-Generate random (probable BPSW) prime of bitsize > 3, pt: prime type of p}

procedure mp_reverse(const a: mp_int; radix: mp_digit; var b: mp_int);
  {-Calculate b with the reversed digits of a in base radix (2..MAXRadix)}

procedure mp_rnr(const a,m: mp_int; var x,y: mp_int);
  {-Rational number reconstruction: for m>0 calculate x,y with a*x=y mod m,}
  { gcd(x,y)=1, 0<=x,|y|<sqrt(m/2), x<>0. x=y=0 if no (unique) solution exists.}

procedure mp_rnr2(const a,m,NN,DD: mp_int; var n,d: mp_int);
  {-Rational number reconstruction: for m,NN,DD > 0 calculate co-prime d,n}
  { with a*d=n mod m, |n|<=N, 0<d<=DD, i.e. a=n/d mod m. Return d=n=0 if }
  { no solution exists. The reconstruction is unique if 2*NN*DD < m.}

function  mp_rqffu(d: longint; var u,v: mp_int): integer;
  {-Return the fundamental unit e = (u + v*sqrt(d))/2 of the real quadratic field}
  { with discriminant d>4, d mod 4 = 0,1. Result is the norm +/- 1 or 0 if error.}

procedure mp_safeprime(var n: mp_int);
  {-Compute the next safe prime >= n: new n and (n-1)/2 are prime, 5 if n<=5}

procedure mp_safeprime_ex(var p,g: mp_int; smallest: boolean);
  {-Compute the next safe prime >= p and a generator g of Z_p^*, smallest g or random}

procedure mp_sigmak(k,n: longint; var a: mp_int);
  {-Compute the sum of the kth powers of the divisors of n; zero if k<0 or n=0.}
  { Special cases k=0: number of divisors of n; k=1: sum of the divisors of n.}

procedure mp_small_factor(const a: mp_int; f0,fmax: mp_digit; var f: mp_digit);
  {-Compute small digit prime factor or 0, f0..fmax, f will be <= min(fmax,$7FFF)}

procedure mp_sqrtmod(const a,p: mp_int; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p, p prime, with Jacobi check.}
  { Err=-1 if (a|p)<>1, Err=1 if failure for p=1 mod 8, Err=2 if p even and <>2}

procedure mp_sqrtmod_ex(const a,p: mp_int; ChkJ: boolean; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p, p prime.}
  { If ChkJ=true then check Jacobi symbol (a|p)=1, Err=-1 if check fails.}
  { Err=1 if failure for p=1 mod 8, Err=2 if p even and not 2.}

procedure mp_sqrtmod2k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate a square root b of an integer a with b*b = a mod 2^k.}
  { Return Err=-1 if there is no solution. For odd a the solutions }
  { are normalized to 0 < b < 2^k/4 for k>3.}

procedure mp_sqrtmodp2(const a,p: mp_int; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p^2, p prime.}
  { Alias for mp_sqrtmodpk(,,2,,). Err: error code from mp_sqrtmodpk.}

procedure mp_sqrtmodpk(const a,p: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate square root b < p^k of a with b*b = a mod p^k, p prime.  }
  { Error codes: if p=2: Err=1 if a<0;  Err=-1 if there is no solution }
  { if p<>2: Err=-1 if (a|p)<>1, Err=1: failure for p=1 mod 8, Err=2 if}
  {          p is even, Err=4: no inverse mod p^(2^i)}

procedure mp_sqrtmodpq(const a,p,q: mp_int; var x,y: mp_int; var Err: integer);
  {-Calculate square roots +x,-x,+y,-y of a mod (pq); p,q primes}
  { If p=q, x is the root from mp_sqrtmodp2 and y is not computed.}
  { For p<>q: Err=4 if gcd(p,q)<>1, otherwise error code from mp_sqrtmod}

function  mp_squad_mod(const a,b,c,p: mp_int; var x1,x2: mp_int): integer;
  {-Solve a^x + b*x + c = 0 mod p, p odd prime; returns number of roots: 0,1,2.}
  { Primality of p is not checked, if 2a has no inverse, b*x + c = 0 is solved.}

procedure mp_swing(n: longint; var a: mp_int);
  {-Calculate Luschny swing number a = swing(n) = n!/((n div 2)!)^2}

function  mp_val(const a: mp_int; r: longint): longint;
  {-Return the valuation of a with respect to r, ie. the largest v with}
  { r^v divides a, v=0 if r=0,1,or -1, v=MaxLongint if a=0}

procedure mp_valrem(const a: mp_int; r: longint; var b: mp_int; var v: longint);
  {-Return valuation v of a with respect to r and remainder b: a = b*r^v}

{#Z+}
{---------------------------------------------------------------------------}
{- 'Internal' functions, don't use them unless you know what you are doing -}
{---------------------------------------------------------------------------}
{#Z-}

function  mp_sqrtmod14_mueller(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 4, return success status.}
  { *internal* assumes a>1,p>1, no argument checks}

function  mp_sqrtmod18_shanks_intern(const a,p: mp_int; var b,q: mp_int; k: longint): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status.}
  { *internal* q and k must be setup with p-1 = 2^k*q}

function  mp_sqrtmod18_shanks(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status; *internal*}

function  mp_sqrtmod18_lucas(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status}
  { *internal* assumes a,p>1, no arg check}

procedure mp_sqrtmod34(const g,p: mp_int; var z: mp_int);
  {-Calculate z with z*z = g mod p, p prime > 0, p mod 4 = 3; *internal*}

procedure mp_sqrtmod58(const g,p: mp_int; var z: mp_int);
  {-Calculate z with z*z = g mod p, p prime = 5 mod 8; *internal*}

function  mp_sqrtmod916_kong(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 9 mod 16, return success status}
  { *internal* assumes a,p>1, no arg check}

procedure s_mp_binom_l(n: longint; k: word; var a: mp_int);
  {-Internal binomial coefficient for small k, no init check}

procedure s_mp_gbinom(n,m,k: longint; var a: mp_int);
  {-Calculate the generalized binomial coefficient a = n!/m!/k!, using prime}
  { power decomposition of a. a must an integer and 0 <= m <= n, 0 <= k <= n}

procedure s_mp_recsplit_fact(n: word; var a: mp_int);
  {-Calculate a = factorial(n) using Recursive Split, error if n > MaxFact}

procedure s_mp_borsch_fact(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) with Borwein/Schoenhage prime factorization method}

procedure s_mp_cbrtmodpk(const a,p: mp_int; k: longint; red: boolean; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod p^k, p prime <> 3, k >= 0.}
  { If red=true, reduce b mod p^k, otherwise b may be >= p^k. Error codes:}
  { 1<Err<4: from mp_cbrtmod, Err=4 if p=3, Err=5: no inverse mod p^(2^i) }

procedure s_mp_cornacchia_ex(const a,b,p,q: mp_int; k: longint; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = m, m=p^k if k>0 otherwise m=p*q, a,b > 0, a+b < m,}
  { (a,m)=1 with Cornacchia's algorithm. Check @x<>@y, but not p,q prime.   }
  { Err <> 0 if error or no solution exist; normalisation x >= y if a=b.    }

function  s_mp_is_cubres(const a, p: mp_int): boolean;
  {-Simple test if a is a cubic residue mod p, p prime. Primality of p}
  { is not checked, but some trivial non-prime cases are detected.}

function  s_mp_is_lpsp(const n: mp_int; strong: boolean): boolean;
  {-Lucas (strong) probable prime test for n using Selfridge A method}

function  s_mp_is_pprime_ex(const a: mp_int; smin,smax: mp_digit): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32); trial division}
  { from smin to smax; no init check, no check if a is less than 2<31 }

function  s_mp_is_pth_power(const a: mp_int; p: longint; var r: mp_int): boolean;
  {-Return true if a is pth power, then a=r^p. a>0, p>2 prime, no checks}

function  s_mp_is_psp2(const n: mp_int): boolean;
  {-Test if n>2 is a base-2 (Fermat) probable prime, no init check}

procedure s_mp_lucasvmod1(const p,n,k,mu: mp_int; var v: mp_int);
  {-Calculate v=v[k] mod n of Lucas V sequence for p,q=1. mu: Barrett parameter}
  { for n, k>0. Internal use: no init check, assumes p<n, return v=2 for k<=0.}

procedure s_mp_npfirst(var a: mp_int; var idx: integer);
  {-Setup idx and increment a to first prime candidate, no init check, a>7}

procedure s_mp_npnext(var a: mp_int; var idx: integer);
  {-Update idx and increment a to next prime candidate, no init check, a>7}

procedure s_mp_nextprime_sieve(safe: boolean; var n: mp_int);
  {-Compute the next prime >= n using a segmented sieve, safe prime if safe=true}

procedure s_mp_nroot_modp(const a,n,p: mp_int; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod p; p prime (not checked), gcd(n, p-1)=1 or gcd(n, (p-1)/n)=1}

procedure s_mp_nroot_modpq(const a,n,p,q: mp_int; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod pq; p,q primes (not checked), gcd(n, phi(pq))=1; gcd(a,p)=1 if p=q}

procedure s_mp_nroot_modpk(const a,n,p: mp_int; k: longint; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod p^k; p prime (not checked), k>0, gcd(a,p)=1 and 1/n mod p^k exists}

procedure s_mp_pell4(const d: mp_int; var x,y: mp_int; var r: integer);
  {-Calculate the smallest non-trivial solution of  x^2 - d*y^2 = r}
  { r in [-4,+4,-1,+1]; returns r=0 if d<2 or d is a square. Uses  }
  { continued fraction expansion of sqrt(d) from Forster [4], ch.25}
  { 'Internal' procedure: Assumes d,x,y are initialized and @x<>@y.}

function  s_mp_rqffu(const d: mp_int; var x,y: mp_int): integer;
  {-Return the fundamental unit e = (u + v*sqrt(d))/2 of the real quadratic field}
  { with discriminant d>4, d mod 4 = 0,1. Result is the norm +/- 1 or 0 if error.}
  { Note: This function uses s_mp_pell4 and is much slower than mp_rqffu but the }
  { discriminant may be large. Beware of very long running times for large d!!   }

procedure s_mp_sqrtmod2k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate unique square root b of an odd integer a with b*b = a mod 2^k and }
  { 0<b<2^k/4 for k<3. Err=1 if a<0; =2 for even a; =3 if a<>1 mod min(2^k,8)}

procedure s_mp_sqrtmodpk(const a,p: mp_int; k: longint; red: boolean; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p^k, p odd prime.}
  { Err=-1 if (a|p)<>1, Err=1 if failure for p=1 mod 8, Err=2 if p even.}
  { Err=4 if no inverse mod p^(2^i). If red=true, reduce b mod p^k, otherwise}
  { b may be >= p^k; ex: a=22,p=3,k=3 --> b=34 and 34^2 = 22 mod 27}

procedure s_mp_valrem_d(const a: mp_int; d: mp_digit; var b: mp_int; var v: longint);
  {-Return valuation v of a with respect to d and remainder b: a = b*d^v}
  { internal iterative version, d no power of two; v=MaxLongint, b=0 if a=0}

implementation


uses
  mp_base, mp_modul, mp_prime;


{---------------------------------------------------------------------------}
procedure s_mp_npfirst(var a: mp_int; var idx: integer);
  {-Setup idx and increment a to first prime candidate, no init check, a>7}
var
  k: integer;
  m: mp_digit;
begin
  idx := -1;
  {make n odd}
  if mp_iseven(a) then mp_inc(a);
  mp_mod_d(a, NPRC_MOD, m);
  if mp_error<>MP_OKAY then exit;
  m := m shr 1;
  {move a to next prime residue class mod MP_MOD and index idx into diff array}
  for k:=m to (NPRC_MOD div 2)-1 do begin
    idx := NPRC_OddIdx[k];
    {note: loop terminates via break because NPRC_OddIdx[(NPRC_MOD div 2)-1]<>-1}
    if idx<>-1 then break;
    mp_add_d(a,2,a);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_npnext(var a: mp_int; var idx: integer);
  {-Update idx and increment a to next prime candidate, no init check, a>7}
begin
  if mp_error<>MP_OKAY then exit;
  {move to next candidate}
  mp_add_d(a,NPRC_Diff[idx],a);
  {get next increment index}
  inc(idx);
  if idx>=NPRC_NRC then idx:=0;
end;


{---------------------------------------------------------------------------}
procedure mp_nextprime(var n: mp_int);
  {-Next prime >= n, 2 if n<=2}
begin
{$ifdef MPC_SieveNextPrime}
  {use sieve for mp_next_prime}
  s_mp_nextprime_sieve(false, n);
{$else}
  mp_nextprime_ex(n, mp_digit(MP_DIGIT_MAX and $3FFF));
{$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_nextprime_ex(var n: mp_int; smax: mp_digit);
  {-Next prime >= n, 2 if n<=2; trial division from 3 to smax}
var
  id: integer;
  l: longint;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_nextprime_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if n<=2 then nextprime = 2}
  if mp_cmp_d(n,3)=MP_LT then begin
    mp_set(n,2);
    exit;
  end;

  {Here n>2. If n <= 2^31-1 then use nextprime32. Note: This is correct}
  {because 2^31-1 is prime and so no 'overflow' to 32 bits can happen!}
  if mp_is_longint(n,l) then begin
    mp_set_int(n,nextprime32(l));
    exit;
  end;

  {move to first candidate}
  s_mp_npfirst(n, id);

  {loop through possible primes}
  repeat
    if mp_error<>MP_OKAY then exit;
    if s_mp_is_pprime_ex(n, 3, smax) then exit;
    {move to next candidate}
    s_mp_npnext(n,id);
  until false;
end;


{---------------------------------------------------------------------------}
procedure mp_prevprime(var n: mp_int);
  {-Previous prime <= n, 0 if n<2}
var
  id,k: integer;
  m: mp_digit;
  l: longint;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_prevprime');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if n<2 then prevprime = 0}
  if mp_cmp_d(n,2)=MP_LT then begin
    mp_zero(n);
    exit;
  end;

  {Here n>=2. If n has less than 32 bits then use prevprime32}
  if mp_is_longint(n,l) then begin
    mp_set_int(n,prevprime32(l));
    exit;
  end;

  {make n odd}
  if mp_iseven(n) then mp_dec(n);

  {$ifndef BIT16}
    {avoid warning, id WILL always be initialized}
    id := 0;
  {$endif}

  mp_mod_d(n, NPRC_MOD, m);
  m := m shr 1;

  {move n to prev prime residue class mod MP_MOD and index id into diff array}
  for k:=m downto 0 do begin
    id := NPRC_OddIdx[k];
    {note: loop is always terminated via break because NPRC_OddIdx[0]<>-1}
    if id<>-1 then break;
    mp_sub_d(n,2,n);
  end;
  repeat
    if mp_error<>MP_OKAY then exit;
    {loop through possible primes}
    if mp_is_pprime(n) then exit;
    {get prev increment index}
    dec(id); if id<0 then id:=NPRC_NRC-1;
    {move to prev candidate}
    mp_sub_d(n,NPRC_Diff[id],n);
  until false;
end;


{---------------------------------------------------------------------------}
procedure s_mp_nextprime_sieve(safe: boolean; var n: mp_int);
  {-Compute the next prime >= n using a segmented sieve, safe prime if safe=true}
const
  MAXSIEVE = 8*4096;
  BitMask  : array[0..7] of byte = ($01,$02,$04,$08,$10,$20,$40,$80);
type
  TBits = array[0..(MAXSIEVE+7) div 8] of byte;
  PBIts = ^TBits;
var
  i,k,p,r,pmax,smax: longint;
  d,f,fmax: mp_digit;
  bsp: PBits;
  ctx: TPrimeContext;
  q: mp_int;
label
  leave;
begin

  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_nextprime_sieve');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n.sign=MP_NEG then mp_zero(n);

  if mp_is_longint(n,i) then begin
    if safe then begin
      i := safeprime32(i);
      if i>0 then begin
        mp_set_int(n,i);
        exit;
      end;
    end
    else begin
      if i<3 then i:=2 else i := nextprime32(i);
      mp_set_int(n,i);
      exit;
    end;
  end;

  if safe then begin
    mp_init(q);
    if mp_error<>MP_OKAY then exit;
  end;

  if mp_iseven(n) then mp_inc(n);

  smax := MAXSIEVE;
  case mp_bitsize(n) of
        0..128: begin
                  pmax := 400;
                  smax := MAXSIEVE div 8;
                end;
      129..256: begin
                  pmax := 800;
                  smax := MAXSIEVE div 4;
                end;
      257..384: begin
                  pmax := 1600;
                  smax := MAXSIEVE div 2;
                end;
      385..512: pmax := 3200;
      513..768: pmax := 6500;
     769..1024: pmax := 15000;
    1025..1280: pmax := 30000;
    1281..1536: pmax := 50000;
    else        pmax := 100000;
  end;
  bsp := mp_alloc(sizeof(TBits));
  if pmax <= longint(mp_max_small) then fmax := pmax
  else fmax := mp_max_small;

  repeat
    fillchar(bsp^, sizeof(TBits), 0);
    findfirstprime32(2,ctx);
    for i:=2 to pmax do begin
      findnextprime32(ctx);
      p := ctx.prime;
      mp_mod_int(n,p,r);
      if mp_error<>MP_OKAY then goto leave;
      if r<>0 then r := p-r;
      while r<SMAX do begin
        k := r shr 3;
        bsp^[k] := bsp^[k] or BitMask[r and 7];
        inc(r,p);
      end;
    end;
    i := 0;
    d := 0;
    while i<SMAX do begin
      if bsp^[i shr 3] and BitMask[i and 7] = 0 then begin
        if d<>0 then mp_add_d(n,d,n);
        d := 0;
        if safe then begin
          mp_sub_d(n,1,q);
          mp_shr1(q);
          {first look for small factors in q}
          mp_small_factor(q, 2, fmax, f);
          if (f=0) then begin
            {then check if n is a psp-2 and q is spsp-2 and slpsp}
            if s_mp_is_psp2(n) then begin
              if mp_is_spsp_d(q,2) and mp_is_slpsp(q) then goto leave;
            end;
          end;
        end
        else begin
          {for n small factors are 'removed' by sieving}
          if mp_is_spsp_d(n,2) and mp_is_slpsp(n) then goto leave;
        end;
      end;
      inc(i,2);
      if d<MP_DIGIT_MAX-2 then inc(d,2)
      else begin
        mp_add_d(n,d+2,n);
        d := 0;
      end;
      if mp_error<>MP_OKAY then goto leave;
    end;
  until false;

leave:
  mp_freemem(pointer(bsp),sizeof(TBits));
  if safe then mp_clear(q);
end;


{---------------------------------------------------------------------------}
procedure mp_safeprime(var n: mp_int);
  {-Compute the next safe prime >= n: new n and (n-1)/2 are prime, 5 if n<=5}
begin
  s_mp_nextprime_sieve(true, n);
end;


{---------------------------------------------------------------------------}
procedure mp_safeprime_ex(var p,g: mp_int; smallest: boolean);
  {-Compute the next safe prime >= p and a generator g of Z_p^*, smallest g or random}
var
  h,t,p1: mp_int;
  d: mp_digit;
const
  MBS = 6;
begin
  if mp_error<>MP_OKAY then exit;
  mp_init3(h,t,p1);
  if mp_error=MP_OKAY then begin
    s_mp_nextprime_sieve(true, p);
    s_mp_sub_d(p,1,p1);
    {h = (p-1) div 2 = p1 div 2 is prime}
    mp_shr(p1,1,h);
    d := 1;
    if mp_bitsize(p) < MBS then smallest := true;
    {Find primitive root mod p. Because the prime factorization of p-1}
    {is 2*h, this reduces to finding 1 < g < p-1 with g^h = -1 mod p. }
    repeat
      if smallest and (d<MP_DIGIT_MAX) then begin
        inc(d);
        if is_square32(d) then inc(d);
        mp_set(g,d);
      end
      else begin
        repeat
          mp_random(p1,g);
        until (mp_error <> MP_OKAY) or (not mp_is_square(g));
      end;
      mp_exptmod(g,h,p,t);
    until (mp_error <> MP_OKAY) or (mp_is_eq(t,p1));
    mp_clear3(h,t,p1);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_lucasvmod1(const p,n,k,mu: mp_int; var v: mp_int);
  {-Calculate v=v[k] mod n of Lucas V sequence for p,q=1. mu: Barrett parameter}
  { for n, k>0. Internal use: no init check, assumes p<n, return v=2 for k<=0.}
var
  i,bc: longint;
  v1: mp_int;
begin

  {return v=2 if k<1}
  bc := mp_bitsize(k);
  if (bc=0) or (k.sign=MP_NEG) then begin
    mp_set(v,2);
    exit;
  end;

  {Initialize v=p, and v1=(p^2-2) mod n if k>1}
  mp_copy(p,v);
  if bc=1 then exit;

  mp_init(v1);            if mp_error<>MP_OKAY then exit;
  mp_sqr(p,v1);           mp_reduce(v1,n,mu);
  mp_sub_d(v1,2,v1);      if v1.sign=MP_NEG then mp_add(v1,n,v1);

  for i:=bc-2 downto 1 do begin
    if mp_isbit(k,i) then begin
      {v = v*v1 - p (mod n), v1 = v1^2 - 2 (mod n)}
      mp_mul(v,v1,v);     mp_reduce(v,n,mu);
      mp_sub(v,p,v);      if v.sign=MP_NEG then mp_add(v,n,v);
      mp_sqr(v1,v1);      mp_reduce(v1,n,mu);
      mp_sub_d(v1,2,v1);  if v1.sign=MP_NEG then mp_add(v1,n,v1);
    end
    else begin
      {v1 = v*v1 - p (mod n), v = v^2 - 2 (mod n)}
      mp_mul(v,v1,v1);    mp_reduce(v1,n,mu);
      mp_sub(v1,p,v1);    if v1.sign=MP_NEG then mp_add(v1,n,v1);
      mp_sqr(v,v);        mp_reduce(v,n,mu);
      mp_sub_d(v,2,v);    if v.sign=MP_NEG then mp_add(v,n,v);
    end;
  end;
  {Calculate final v (v1 not needed)}
  if mp_isodd(k) then begin
    mp_mul(v,v1,v);       mp_reduce(v,n,mu);
    mp_sub(v,p,v);        if v.sign=MP_NEG then mp_add(v,n,v);
  end
  else begin
    mp_sqr(v,v);          mp_reduce(v,n,mu);
    mp_sub_d(v,2,v);      if v.sign=MP_NEG then mp_add(v,n,v);
  end;
  mp_clear(v1);
end;


{---------------------------------------------------------------------------}
procedure mp_4sq(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}
var
  n,p,x: mp_int;
  n8,i, Err: integer;
  v,s: longint;
label
  done;
const
  sn: array[0..17] of word = (1,2,3,10,34,58,85,130,214,226,370,526,706,730,1414,1906,2986,9634);
  sa: array[0..17] of byte = (1,1,1, 1, 3, 3, 6,  3,  3,  8,  8,  6, 15,  1,   6,  13,  21,  56);
  sb: array[0..17] of byte = (0,1,1, 3, 3, 7, 7, 11,  6,  9,  9,  7, 15, 27,  17,  21,  32,  57);
  sc: array[0..17] of byte = (0,0,1, 0, 4, 0, 0,  0, 13,  9, 15, 21, 16,  0,  33,  36,  39,  57);
begin

  {Based on an algorithm by Peter Schorn http://www.schorn.ch/lagrange.html}
  {and his other resources howto.html, FourSquares.java, and decompose.py}

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(m) or mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_4sq');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {must be positive}
  if m.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_4sq: m < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if mp_is0(m) then begin
    mp_zero(a);
    mp_zero(b);
    mp_zero(c);
    mp_zero(d);
    exit;
  end;

  mp_init3(n,p,x);
  if mp_error<>MP_OKAY then exit;

  {Divide out powers of 4: m = 4^v*n}
  v := mp_cnt_lsb(m) div 2;
  mp_shr(m, 2*v, n);

  {Initialize default values after computing v,n: @m may be in [@a,@b,@c,@d]}
  mp_zero(a);
  mp_zero(b);
  mp_zero(c);
  mp_zero(d);

  {x = sqrt(n), p = n-x^2}
  mp_sqrtrem(n,x,p);
  if mp_is0(p) then begin
    {n is a square, almost done}
    mp_exch(a,x);
    goto done;
  end;

  {get n mod 8}
  n8 := n.pdigits^[0] and 7;
  if n8=7 then begin
    {If n=7 mod 8, then all four variables a,b,c,d<>0 are needed:}
    {Subtract largest square d^2 < n with d^2 mod 8 <> 0}
    if (x.used>0) and (x.pdigits^[0] and 3 <> 0) then begin
      {d = x = floor(sqrt(n)) already found}
      mp_exch(x,d);
      mp_exch(n,p);
    end
    else begin
      {x=0 mod 8, set d=x-1. Then n-d^2 = n-x^2 + 2x-1 = p + 2x-1}
      mp_sub_d(x,1,d);
      mp_shl(x,1,n);
      mp_dec(n);
      mp_add(n,p,n);
    end;
    {multiply d by 4^v, this scaling is skipped for d at label done}
    mp_shl(d,v,d);
    {Recalculate x,p with n = x^2 + p}
    mp_sqrtrem(n,x,p);
  end
  else if n8 and 3 = 1 then begin
    {Here n mod 4=1; if n is prime then n = a^2 + b^2 via Cornacchia}
    if mp_is_pprime(n) then begin
      mp_set(c,1);
      mp_cornacchia(c,n,a,b,Err);
      if Err=0 then begin
        mp_zero(c);
        goto done;
      end;
    end;
  end;

  {Check special cases which cannot be handled by rest of algorithm}
  if mp_is_longint(n,s) and (s<=9634) then begin
    for i:=0 to 17 do begin
      if s=sn[i] then begin
        mp_set(a,sa[i]);
        mp_set(b,sb[i]);
        mp_set(c,sc[i]);
        goto done;
      end;
    end;
  end;

  {Here n = x^2 + p with n<>0, n mod 8 <> 7, n mod 4 <> 0}
  {$ifdef MPC_USE_Assert}
    assert(n.used>0, MPAF+'n<>0');
    assert(n.pdigits^[0] and 7 <> 7, MPAF+'n mod 8 <> 7');
    assert(n.pdigits^[0] and 3 <> 0, MPAF+'n mod 4 <> 0');
  {$endif}

  if 3 = n.pdigits^[0] and 3 then begin
    {Case 3 = n mod 4 = n mod 8; find odd x with p = (n-x^2)/2 prime}
    if mp_iseven(x) then begin
      mp_add(p,x,p);
      mp_dec(x);
      mp_add(p,x,p);
    end;
    mp_shr(p,1,p);
    repeat
      if mp_is_pprime(p) then begin
        mp_set(n,1);
        mp_cornacchia(n,p,b,n,Err);
        if Err=0 then begin
          mp_sub(b,n,c);
          mp_abs(c,c);
          mp_add(b,n,b);
          mp_copy(x,a);
          goto done;
        end;
      end;
      {Next pair (p,x): p = p + 2(x-1);  x = x - 2}
      mp_sub_d(x,1,n);
      mp_shl(n,1,n);
      mp_add(p,n,p);
      mp_sub_d(x,2,x);
    until (mp_error<>MP_OKAY) or (x.sign=MP_NEG);
  end
  else begin
    {Case n mod 4 = 1 or 2; find even x with p = n - x^2 prime}
    if (n.pdigits^[0] and 1)=(x.pdigits^[0] and 1) then begin
      mp_shl(x,1,n);
      mp_add(p,n,p);
      mp_dec(x);
      mp_dec(p);
    end;
    repeat
      if mp_is_pprime(p) then begin
        mp_set(n,1);
        mp_cornacchia(n,p,b,c,Err);
        if Err=0 then begin
          mp_copy(x,a);
          goto done;
        end;
      end;
      {Next pair (p,x): p = p + 4(x-1);  x = x - 2}
      mp_sub_d(x,1,n);
      mp_shl(n,2,n);
      mp_add(p,n,p);
      mp_sub_d(x,2,x);
    until (mp_error<>MP_OKAY) or (x.sign=MP_NEG);
  end;

  {no solution found, reset d=0. a,b,c are still zero}
  mp_zero(d);

done:
  {rescale initial powers of 4, d is 0 or already scaled}
  if v>0 then begin
    mp_shl(a,v,a);
    mp_shl(b,v,b);
    mp_shl(c,v,c);
  end;

  mp_clear3(n,p,x);
end;


{---------------------------------------------------------------------------}
procedure mp_4sq_sa(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2, a<=b<=c<=d.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}
begin
  {Error checking in mp_4sq}
  mp_4sq(m,a,b,c,d);
  { make a <= b <= c <= d}
  if mp_error<>MP_OKAY then exit;
  if mp_cmp_mag(a,b)>0 then mp_exch(a,b);
  if mp_cmp_mag(a,c)>0 then mp_exch(a,c);
  if mp_cmp_mag(a,d)>0 then mp_exch(a,d);
  if mp_cmp_mag(b,c)>0 then mp_exch(b,c);
  if mp_cmp_mag(b,d)>0 then mp_exch(b,d);
  if mp_cmp_mag(c,d)>0 then mp_exch(c,d);
end;


{---------------------------------------------------------------------------}
procedure mp_4sq_sd(const m: mp_int; var a,b,c,d: mp_int);
  {-Decompose m>=0 into 4 squares: m = a^2 + b^2 + c^2 + d^2, a>=b>=c>=d.}
  { If no solutions for m>0 was found, a=b=c=d=0 is returned.}
begin
  {Error checking in mp_4sq}
  mp_4sq(m,a,b,c,d);
  { make a <= b <= c <= d}
  if mp_error<>MP_OKAY then exit;
  if mp_cmp_mag(a,b)<0 then mp_exch(a,b);
  if mp_cmp_mag(a,c)<0 then mp_exch(a,c);
  if mp_cmp_mag(a,d)<0 then mp_exch(a,d);
  if mp_cmp_mag(b,c)<0 then mp_exch(b,c);
  if mp_cmp_mag(b,d)<0 then mp_exch(b,d);
  if mp_cmp_mag(c,d)<0 then mp_exch(c,d);
end;


{---------------------------------------------------------------------------}
procedure s_mp_binom_l(n: longint; k: word; var a: mp_int);
  {-Internal binomial coefficient for small k, no init check}
var
  i,e2,q: word;
  p: longint;
  tmp: mp_int;
begin
  if mp_error<>MP_OKAY then exit;

  {special cases}
  if k>n then begin
    mp_zero(a);
    exit;
  end;
  if k=1 then begin
    mp_set_int(a,n);
    exit;
  end;
  mp_set1(a);
  if (n<2) or (k=0) or (k=n) then exit;

  mp_init(tmp);
  if mp_error<>MP_OKAY then exit;

  {if k>n/2 use symmetry to reduce k}
  if k>n-k then k := n-k;
  n := n-k+1;
  {e2 will count powers of 2}
  e2 := 0;
  for i:=1 to k do begin
    {divide out powers of 2 in numerator}
    p := n;
    while p and 1 = 0 do begin
      p := p shr 1;
      inc(e2);
    end;
    s_mp_mul_int(a,p,a,tmp);
    {divide out powers of 2 in denominator}
    q := i;
    while q and 1 = 0 do begin
      q := q shr 1;
      dec(e2);
    end;
    mp_div_w(a,q,@a,q);
    inc(n);
  end;
  if e2>0 then mp_shl(a,e2,a);
  mp_clear(tmp);
end;


{---------------------------------------------------------------------------}
function GBinPrimeExpo(p,n,m,k: longint): longint;
  {-Compute max. exponent e with p^e divides the integer n!/m!/k!, n >= m >= k}
var
  e: longint;
begin
  e := 0;
  while k >= p do begin
    n := n div p;
    m := m div p;
    k := k div p;
    inc(e, n - m - k);
  end;
  while m >= p do begin
    n := n div p;
    m := m div p;
    inc(e, n - m);
  end;
  while n >= p do begin
    n := n div p;
    inc(e, n);
  end;
  GBinPrimeExpo := e;
end;


{---------------------------------------------------------------------------}
procedure s_mp_gbinom(n,m,k: longint; var a: mp_int);
  {-Calculate the generalized binomial coefficient a = n!/m!/k!, using prime}
  { power decomposition of a. a must an integer and 0 <= m <= n, 0 <= k <= n}

  {The basic idea is from Peter Luschny's [11] FastBinomial page }
  {http://www.luschny.de/math/factorial/FastBinomialFunction.html}

  {For each prime 2 < p <= n the power p^e dividing  n!/m!/k! is }
  {computed, see P.Goetgheluck, Computing Binomial Coefficients, }
  {Amer. Math. Monthly 94 (1987), 360-365. See also H. Reddmann's}
  {public domain source NCombi.pas in DEC 5, Part II, available  }
  {from http://code.google.com/p/delphidec/                      }

  {Also used is Marcel Martin's [8] idea to evaluate the product }
  {of a list of numbers: It is faster to calculate intermediate  }
  {products of roughly the same size rather sequentially multiply}
  {each factor of the list. The arrays X and Level are used for  }
  {managing the intermediates with a stack like structure.       }

const
  MAXS = 31;

var
  X    : array[0..MAXS] of mp_int;   {Used for managing the intermediate}
  Level: array[0..MAXS] of integer;  {products of the prime powers p^e}

var
  e,p,e2: longint;
  s,smax: integer;
  {$ifdef BIT16}
    ctx: TPrimeContext;
  {$else}
    sieve: TSieve;
  {$endif}

label
  leave;

begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_gbinom');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if (m<0) or (k<0) or (m>n) or (k>n) then begin
    {$ifdef MPC_HaltOnArgCheck}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mp_gbinom: invalid arguments');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {init default result}
  mp_set1(a);

  {Init partial product stack indices and some derived variables}
  s   := -1;
  smax:= -1;

  {force m >= k}
  if k>m then begin
    e := k;
    k := m;
    m := e;
  end;

  {powers of 2 are calculated separately and shifted in at the end}
  e2 := (n-popcount32(n)) - (k-popcount32(k)) - (m-popcount32(m));
  if e2<0 then begin
    {$ifdef MPC_HaltOnArgCheck}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mp_gbinom: negative exponent for p=2');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {$ifdef BIT16}
    {Init context so that FindNextPrime32 returns 3}
    FindFirstPrime32(2,ctx);
  {$else}
    prime_sieve_init(sieve, 3);
  {$endif}

  {iteration over the primes > 2}
  while mp_error=MP_OKAY do begin
    {$ifdef BIT16}
      FindNextPrime32(ctx); p := ctx.prime;
    {$else}
      p := prime_sieve_next(sieve);
    {$endif}

    {done if all primes <= n are processed}
    if p>n then break;

    {get max exponent for e with p^e | n!/m!/k!}
    e := GBinPrimeExpo(p,n,m,k);
    if e<0 then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('s_mp_gbinom: negative prime exponent');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;

    {if p divides n!/m!/k!, accumulate p^e}
    if e>0 then begin
      if s<MAXS then begin
        inc(s);
        if s>smax then begin
          {init a new stack level, smax is only incremented if no error}
          mp_init(X[s]);
          if mp_error<>MP_OKAY then goto leave;
          smax := s;
        end;
        {calculate and push p^e }
        if e=1 then mp_set_int(X[s],p)
        else if (e=2) and (p<46341) then mp_set_int(X[s],sqr(p))
        else mp_set_pow(X[s],p,e);
        {new factor has level 0}
        Level[s] := 0;
        {reduce stack}
        while (s > 0) and (Level[s-1] = Level[s]) do begin
          mp_mul(X[s-1],X[s],X[s-1]);
          dec(s);
          {inc level of new partial factor}
          inc(Level[s]);
        end;
      end
      else begin
       {MAXS should be big enough, just paranoia}
       {$ifdef MPC_HaltOnError}
         {$ifdef MPC_UseExceptions}
           raise MPXRange.Create('s_mp_gbinom: stack too small');
         {$else}
           RunError(MP_RTE_RANGE);
         {$endif}
       {$else}
         set_mp_error(MP_RANGE);
         goto leave;
       {$endif}
      end;
    end;

  end;
  if (mp_error=MP_OKAY) and (s>=0) then begin
    {finalize partial product processing if stack is not empty}
    mp_exch(a,X[s]);
    while s > 0 do begin
      dec(s);
      mp_mul(a,X[s],a);
    end;
  end;
  {shift in powers of 2}
  mp_shl(a,e2,a);

leave:
  {$ifndef BIT16}
    prime_sieve_clear(sieve);
  {$endif}
  {clear only the used X[s]}
  for s:=0 to smax do mp_clear(X[s]);
end;


{---------------------------------------------------------------------------}
procedure mp_binomial(n,k: longint; var a: mp_int);
  {-Calculate the binomial coefficient a = (n choose k)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_binomial');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {init default result}
  mp_set1(a);
  if (k=0) or (k=n) then exit;

  if (k=1) or (k=pred(n)) then begin
    mp_set_int(a,n);
    exit;
  end;

  if k<0 then begin
    if (n>=0) or (n<k) then mp_zero(a)
    else begin
      {Here k <= n < 0:  use identity 1.2.6 (19) from Knuth[3]}
      {binomial(n,k) = (-1)^(n-k)*binomial(n-k-n-1,n-k)}
      n := n-k;
      mp_binomial(-succ(k),n,a);
      if odd(n) then s_mp_chs(a);
    end;
    exit;
  end;

  if n<0 then begin
    {Here n < 0 and k >= 0, use Knuth[3], 1.2.6 (17)}
    {binomial(n,k) = (-1)^k * binomial(k-n-1,k)}
    mp_binomial(k-n-1,k,a);
    if odd(k) then s_mp_chs(a);
    exit;
  end;

  {here n >= 0, handle easy case first}
  if k>n then begin
    mp_zero(a);
    exit;
  end;

  {Here n >= k >= 2. If k>n/2 use symmetry to reduce k}
  if k>n-k then k := n-k;

  if k<64 then begin
    {handle small values without prime power method}
    s_mp_binom_l(n,word(k),a)
  end
  else begin
    {Use generalized method with prime powers}
    s_mp_gbinom(n,n-k,k,a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_catalan(n: longint; var c: mp_int);
  {-Return the n'th Catalan number C_n = binomial(2n,n)/(n+1), 0 if n < 0}
var
  d: mp_digit;
begin
  {MaxCatalan ~ MP_Maxbit/2.0005}
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_catalan');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n<0 then mp_zero(c)
  else if n<64 then begin
    {C_n = binomial(2n,n)/(n+1)}
    s_mp_binom_l(n+n,n,c);
    mp_div_d(c,n+1,@c,d);
  end
  else begin
    {Division by n+1 would be very slow for MP_16BIT and n >= 2^15, because}
    {it will be a multiprecision division, use the following form for C_n: }
    {C_n = (2n)!/n!/(n+1)!}
    s_mp_gbinom(n+n,n+1,n,c);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_cbrtmod_ex(const a, p: mp_int; var x: mp_int; pru: pmp_int; var Err: integer);
  {-Compute a cube root x of a with x^3 = a mod p, p prime. If pru<>nil, pru^ }
  { is set to a 3rd root of unity (or 1 if x is unique), the other roots are  }
  { x*pru^, x*pru^*pru^. Err <> 0, if a is not a 3rd power or p is not prime. }
var
  i,j,k: longint;
  d: mp_digit;
  b,q,t,y,z: mp_int;
label
  leave;
begin

  {See S.C. Lindhurst [32], Chap. 3.2: Extension of Shanks's Algorithm to rth Roots}
  Err := 0;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(x) or ((pru<>nil) and mp_not_init(pru^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cbrtmod_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (@a=@p) or (@a=@x) or (@a=pru) or (@p=@x) or (@p=pru) or (@x=pru) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_cbrtmod_ex: var addresses not unique');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  if p.used=1 then begin
    d := p.pdigits^[0];
    if (d=2) or (d=3) then begin
      if pru<>nil then mp_set1(pru^);
      mp_mod_d(a,d,d);
      mp_set(x,d);
      exit;
    end;
  end;

  if s_mp_is_le0(p) or mp_is1(p) or mp_iseven(p) then begin
    {p<2 or even, ie p is not prime}
    Err := 1;
    exit;
  end;

  mp_init5(b,q,t,y,z);
  if mp_error<>MP_OKAY then exit;

  {case a=0 mod p}
  mp_mod(a,p,t);
  if t.used=0 then begin
    mp_zero(x);
    if pru<>nil then mp_set1(pru^);
    goto leave;
  end;

  mp_mod_d(p,3,d);

  if d=0 then begin
    {p>3 is a multiple of 3, i.e. p is not prime}
    Err := 1;
    goto leave;
  end
  else if d=2 then begin
    {Easy case: if p=2 mod 3, then  x = a^((2p-1)/3) mod p}
    mp_add(p,p,t);
    mp_dec(t);
    mp_div_d(t,3,@t,d);
    mp_exptmod(a,t,p,x);
    if pru<>nil then mp_set1(pru^);
    goto leave;
  end;

  {Here p=1 mod 3; calculate n,q with p-1 = 3^k*q and gcd(q,3)=1}

  {Compute b=(p-1)/3: used for k,q calculation, non-3rd power search,}
  {primitive 3rd root determination, and main loop initialization}
  mp_sub_d(p,1,b);
  mp_div_d(b,3,@t,d);

  {Get k and q, because p=1 mod 3 we can use repeat}
  k := 0;
  mp_copy(t,b);
  repeat
    inc(k);
    mp_exch(q,t);
    mp_div_d(q,3,@t,d)
  until d<>0;

  {find a cubic non-residue z}
  i := 1;
  repeat
    i := nextprime32(i+1);
    if i<0 then begin
      Err := 2;
      goto leave;
    end;
    mp_set_int(z,i);
    mp_exptmod(z,b,p,y);
    if not mp_is1(y) then break;
    if mp_error<>MP_OKAY then goto leave;
    if i=541 then begin
      {after unsuccessfully testing 100 primes as CNR candidate,}
      {check if p itself is prime before going on}
      if not mp_is_pprime(p) then begin
        Err := 1;
        goto leave;
      end;
    end;
  until false;

  {z is a cubic non-residue,  y = z^((p-1)/3) mod p is primitive 3rd root}
  {of unity. Initialize loop: (k=n), z = z^q, x = a^((mq+1)/3), b = a^mq,}
  {where m=1 if q mod 3=2, and m=2 if q mod 3 = 1}

  {z = z^q}
  mp_exptmod(z,q,p,z);

  {t = (mq-2)/3}
  mp_mod_d(q,3,d);
  if d=1 then mp_shl1(q);
  mp_sub_d(q,2,t);
  mp_div_d(t,3,@t,d);

  {x = a^((mq-2)/3)}
  mp_exptmod(a,t,p,x);

  {b = x*(ax)^2 = a^2*x^3 = a^(2+3(mq-2)/3) = a^mq}
  mp_mulmod(a,x,p,b);
  mp_sqrmod(b,p,b);
  mp_mulmod(b,x,p,b);

  {x = x*a = a^(1+(mq-2)/3) = a^((mq+1)/3)}
  mp_mulmod(a,x,p,x);

  {main loop, invariant: x^3 = a*b}
  while not mp_is1(b) do begin
    if mp_error<>MP_OKAY then goto leave;
    {find least positive integer with b^(3^j) = 1}
    j := 0;
    mp_copy(b,q);
    repeat
      inc(j);
      {j=k should not happen if p is prime and a is a 3rd power, indicate failure}
      if (j=k) or (mp_error<>MP_OKAY) then begin
        Err := 3;
        goto leave;
      end;
      mp_sqrmod(q,p,t);
      mp_mulmod(q,t,p,t);
      if mp_is1(t) then break;
      mp_exch(q,t);
    until false;

    {Get l=1 or 2 with (bz^l)^(3^j-1) = b^(3^j-1)*(z^(3^j-1))^l = 1 }
    {c.f.3.2.3 Speeding Things Up. We have q=b^(3^j-1), q^3=1, q<>1,}
    {i.e. q=y or q=y^2. If q=y, then l=2; continue with z^2 and y^2.}
    {see also cubic_root.py from http://tnt.math.se.tmu.ac.jp/nzmath}
    if mp_is_eq(q,y) then begin
      mp_sqrmod(z,p,t);
      mp_sqrmod(y,p,y);
    end
    else mp_copy(z,t);

    {t = t^(3^(k-j-1))}
    mp_set(q,3);
    mp_exptmod_d(q,k-j-1,p,q);
    mp_exptmod(t,q,p,t);

    {z = t^3}
    mp_sqrmod(t,p,z);
    mp_mulmod(z,t,p,z);
    {x = x*t}
    mp_mulmod(x,t,p,x);
    {b = b*z}
    mp_mulmod(b,z,p,b);
    k := j;
  end;

  if pru<>nil then mp_copy(y,pru^);

leave:
  mp_clear5(b,q,t,y,z);
end;


{---------------------------------------------------------------------------}
procedure mp_cbrtmod(const a, p: mp_int; var x: mp_int; var Err: integer);
  {-Compute a cube root x of a with x^3 = a mod p, p prime.  }
  { Err in [1..3], if a is not a 3rd power or p is not prime.}
begin
  mp_cbrtmod_ex(a, p, x, nil, Err);
end;


{---------------------------------------------------------------------------}
procedure s_mp_cbrtmodpk(const a,p: mp_int; k: longint; red: boolean; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod p^k, p prime <> 3, k >= 0.}
  { If red=true, reduce b mod p^k, otherwise b may be >= p^k. Error codes:}
  { 1<Err<4: from mp_cbrtmod, Err=4 if p=3, Err=5: no inverse mod p^(2^i) }
var
  p2i,ri,x,z: mp_int;
  i,im: integer;
begin
  if k<0 then begin
    Err := 99;
    exit;
  end;
  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_cbrtmodpk');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_iszero(a) or mp_is1(a) or (k=0) then begin
    {if a=0 or 1, set b=a}
    if k=0 then mp_zero(b) else mp_copy(a,b);
    exit;
  end;
  if mp_cmp_d(p,3)<>0 then begin
    mp_init4(p2i,ri,x,z);
    if mp_error<>MP_OKAY then exit;
    {r0 = cbrt(a) mod p}
    mp_cbrtmod(a,p,ri,Err);
    if (Err=0) and (k>1) then begin
      if mp_is0(ri) and (k>=3) then begin
        {zero solution but a<>0: try to divide out p^3}
        mp_sqr(p,p2i);
        mp_mul(p2i,p,p2i);
        mp_divrem(a,p2i,@x,@z);
        if mp_is0(z) then begin
          s_mp_cbrtmodpk(x,p,k-3,true,ri,Err);
          if Err=0 then mp_mul(ri,p,ri);
        end
        else Err:=5;
      end
      else begin
        {[10] Alg 2.3.11 (Hensel lifting) with f(x)=a-x^3, f'(x)=-3x}
        {Calls function newr repeatedly until 2^i >= k; here k>1    }
        im := bitsize32(k)-1;
        {2^im <= k <= 2^(im+1)}
        if k and pred(k)=0 then begin
          {k a power of two: decrement upper bound}
          dec(im);
          {Now k=2^(im+1) and we are done after im+1 loops}
        end;
        for i:=0 to im do begin
          {calculate p^(2^i): copy from p or square p^(2^(i-1))}
          if i=0 then mp_copy(p,p2i) else mp_sqr(p2i,p2i);
          {z = f'(ri)^-1 mod p^(2^i) = -(3*ri^2)^-1 mod p^(2^i)}
          mp_sqr(ri,z);
          mp_mul_d(z,3,z);
          if not mp_invmodf(z,p2i,z) then begin
            Err := 8;
            break;
          end;
          {x = f(ri)/p^(2^i) = (a-ri^3)/p^(2^i)}
          mp_sqr(ri,x);
          mp_mul(x,ri,x);
          mp_sub(a,x,x);
          mp_div(x,p2i,x);
          {x = -xz mod p^(2^i)) = x/(3*ri^2) mod p^(2^i))}
          mp_mulmod(z,x,p2i,x);
          {r(i+1) = ri + x*p^(2^i)}
          mp_mul(x,p2i,x);
          mp_add(x,ri,ri);
        end;
      end;
    end;
    if Err=0 then begin
      {store last ri as result into b, reduce if requested}
      if red then begin
        mp_expt_int(p,k,x);
        mp_mod(ri,x,b);
      end
      else mp_exch(b,ri);
    end;
    mp_clear4(p2i,ri,x,z);
  end
  else Err:=4;
end;


{---------------------------------------------------------------------------}
function s_mp_is_cubres(const a, p: mp_int): boolean;
  {-Simple test if a is a cubic residue mod p, p prime. Primality of p}
  { is not checked, but some trivial non-prime cases are detected.}
var
  d: mp_digit;
  t: mp_int;
begin
  s_mp_is_cubres := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_is_cubres');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if p.sign=MP_NEG then exit;

  if p.used=1 then begin
    d := p.pdigits^[0];
    if d<2 then exit;  {p=0 or 1}
    if d<4 then begin
      {p=2 or 3: all a are CR}
      s_mp_is_cubres := true;
      exit;
    end;
  end;
  if mp_iseven(p) then exit;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_mod(a,p,t);
  if mp_is0(t) or mp_is1(t) then begin
    {case a=0,1 mod p}
    s_mp_is_cubres := true
  end
  else begin
    {Check t=p-1=-1 mod p. Don't waste time for inc(t) if t is too small.}
    {Here 0 <= t < p; and since p is odd, t.used <= p.used = (p-1).used}
    if t.used >= p.used then begin
      mp_inc(t);
      if mp_cmp_mag(p,t)=0 then begin
        mp_clear(t);
        s_mp_is_cubres := true;
        exit;
      end;
    end;
    mp_sub_d(p,1,t);
    mp_div_d(t,3,@t,d);
    if d=0 then begin
      {Euler's criterion for p=1 mod 3: x^3 = a mod is solvable if }
      {a^((p-1)/gcd(p-1,3)) = a^((p-1)/3)) = a^t = 1 mod p}
      mp_exptmod(a,t,p,t);
      s_mp_is_cubres := mp_is1(t);
    end
    else begin
      {Easy cases: if d=2 ie p=0 mod 3, then p is not prime and  }
      {if d=1 ie p=2 mod 3, then a^((2p-1)/3) mod p is a solution}
      s_mp_is_cubres := d=1;
    end;
  end;
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_cbrtmod3k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod 3^k. Error codes: 1<Err<4 code from}
  { mp_cbrtmod, Err=4 if p=3, Err=5: no inverse mod 3^(2^i), Err=6: (internal d<>0)}
var
  i: integer;
  d: mp_digit;
  p2i,pk,x,z: mp_int;
begin
  if k<=0 then begin
    Err := 99;
    exit;
  end;

  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cbrtmod3k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_iszero(a) or mp_is1(a) then begin
    {if a=0 or 1, set b=a}
    mp_copy(a,b);
    exit;
  end;
  mp_init4(p2i,pk,x,z);
  if mp_error<>MP_OKAY then exit;
  mp_set_pow(pk,3,k);  {3^k}
  mp_set(p2i,3);       {3^2^i}
  {get solution mod 3}
  mp_cbrtmod(a,p2i,z,Err);
  {if no solution mod 3 then no solution mod 3^k}
  if (Err=0) and (k>1) then begin
    if mp_is0(z) and (k>=3) then begin
      {zero solution but a<>0: try to divide out 3^3}
      mp_div_d(a,27,@x,d);
      if d=0 then begin
        mp_cbrtmod3k(x,k-3,z,Err);
        if Err=0 then mp_mul_d(z,3,z);
      end
      else Err:=5;
    end
    else begin
      {loop is executed until 3^(2^i) >= 3^k, i.e. 2^i>=k}
      for i:=0 to bitsize32(k)-1 do begin
        {z<>0 is a solution mod 3^(2^i), calculate solution mod 3^(2^(i+1))}
        {s_mp_cbrtmodpk cannot be used because f'(z)= 0 mod 3. Use}
        {discrete Newton z~ = z - f(z)/f'(z)   = z - (z^3-a)/3z^2 }
        {                   = z - z/3 + a/3z^2 = (2z + a/z^2)/3   }
        {i.e. iterate z = (a*z^-2 +2z)/3 mod 3^(2^(i+1))}
        mp_sqr(p2i,p2i);
        mp_sqr(z,x);
        if not mp_invmodf(x,p2i,x) then begin
          Err := 5;
          break;
        end;
        mp_mul(a,x,x);
        mp_shl1(z);
        mp_add(z,x,z);
        mp_div_d(z,3,@z,d);
        if d<>0 then begin
          Err := 6;
          break;
        end;
        mp_mod(z,p2i,z);
      end;
    end;
  end;
  if Err=0 then begin
    {z is a solution, either directly for k=1 or via iteration for k>1}
    {always reduce mod 3^k}
    mp_mod(z,pk,b);
  end;
  mp_clear4(p2i,pk,x,z);
end;


{---------------------------------------------------------------------------}
procedure mp_cbrtmodpk(const a,p: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Compute a cube root b of a with b^3 = a mod p^k, p prime. b will be  }
  { reduced mod p^k. Error codes: 1<Err<4: from mp_cbrtmod, Err=4 if p=3,}
  { Err=5: no inverse mod p^(2^i), Err=6: internal check for p=3 failed. }
begin
  if mp_cmp_d(p,3)=0 then mp_cbrtmod3k(a,k,b,Err)
  else s_mp_cbrtmodpk(a,p,k,true,b,Err);
end;


{---------------------------------------------------------------------------}
procedure mp_cbrtmodpq(const a,p,q: mp_int; var x: mp_int; var Err: integer);
  {-Compute a cube root x of a mod (pq); p,q primes. If p=q, mp_cbrtmodpk is}
  { used. For p<>q: Err=8 if gcd(p,q)<>1, otherwise error code from mp_cbrtmod.}
var
  c,d,n,r: mp_int;
begin
  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(q) or mp_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cbrtmodpq');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Handle case p=q}
  if mp_is_eq(p,q) then begin
    mp_cbrtmodpk(a,p,2,x,Err);
    exit;
  end;

  mp_init4(c,d,n,r);
  {First compute extended gcd, c*p + d*q = gcd(p,q) = n. Error if n<>1}
  {r = cbrt(a) mod p, s = cbrt(a) mod q. Result x = (rdq+scp) mod pq}
  mp_xgcd(p,q,@c,@d,@n);
  if not mp_is1(n) then Err:=8
  else begin
    {calculate roots mod p and mod q}
    mp_cbrtmod(a,p,r,Err);
    if (Err=0) and (mp_error=MP_OKAY) then begin
      mp_mul(p,q,n);
      mp_mulmod(d,q,n,d);
      mp_mulmod(r,d,n,d);
      mp_cbrtmod(a,q,r,Err);
      if (Err=0) and (mp_error=MP_OKAY) then begin
        mp_mulmod(c,p,n,c);
        mp_mulmod(r,c,n,c);
        mp_addmod(c,d,n,x);
      end;
    end;
  end;
  mp_clear4(c,d,n,r);
end;


{---------------------------------------------------------------------------}
procedure mp_cornacchia(const d,p: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve x^2 + d*y^2 = p, p prime, 0<d<p with Cornacchia's algorithm. Check}
  { @x<>@y, but not p prime. Err<>0 if no solution exist. x >= y if d=1.}
var
  a,b,t: mp_int;
  d1: boolean;
begin
  {Ref: [24] Cohen, Algorithm 1.5.2}
  {     [10] Crandall/Pomerance, Algorithm 2.3.12}

  {Error codes:
    -2:  d<=0 or d>p or initial mp_error<>0
     4:  (p-x^2)/d is no square
     5:  (p-x^2)/d is no integer
   else  error code from mp_sqrtmod}

  Err := -2;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(d) or mp_not_init(p) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cornacchia');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @x=@y then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_cornacchia: @x=@y');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  if s_mp_is_le0(d) then exit;

  mp_init3(a,b,t);

  {Check upper bound for d and handle special cases}
  mp_add_d(d,1,t);
  if mp_is_ge(t,p) then begin
    {no solution for d>p}
    if mp_is_le(d,p) then begin
      {d=p or d=p-1}
      Err := 0;
      mp_set1(y);
      if mp_is_eq(d,p) then begin
        {d=p  => x=0, y=1}
        mp_zero(x);
      end
      else begin
        {d=p-1  => x=1, y=1}
        mp_set1(x);
      end;
    end;
    mp_clear3(a,b,t);
    exit;
  end;

  {Remember if d=1}
  d1 := mp_is1(d);

  {Step 1a: Solve b^2 = -d mod p, Jacobi check (-d|p) is done in mp_sqrtmod}
  mp_chs(d,t);
  mp_sqrtmod(t,p,b,Err);
  if Err=0 then begin
    {Step 1b: if b <= p/2, or equivalent 2b <= p, then replace b by p-b}
    mp_mul_2(b,t);
    if mp_is_le(t,p) then mp_sub(p,b,b);
    mp_copy(p,a);
    mp_sqrt(p,t);
    {Step 2: Apply Euclidean algorithm to (a,b) until b <= sqrt(p)}
    while (mp_error=MP_OKAY) and mp_is_gt(b,t) do begin
      mp_mod(a,b,a);
      mp_exch(a,b);
    end;
    {Step 3: Set x=b, y = sqrt((p-x^2)/d). If y is in N, then the}
    {solution is (x,y). Otherwise there is no solution, return Err}
    mp_sqr(b,t);
    mp_sub(p,t,a);
    mp_divrem(a,d,@a,@t);
    if mp_is0(t) then begin
      if mp_is_square2(a,@y) then begin
        mp_exch(x,b);
        {Normalize solution: x >= y if d=1}
        if d1 and mp_is_lt(x,y) then mp_exch(x,y);
      end
      else Err := 4;
    end
    else Err := 5;
  end;

  mp_clear3(a,b,t);
end;


{---------------------------------------------------------------------------}
procedure mp_cornacchia4(const d,p: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve x^2 + |d|*y^2 = 4p, p prime, -4p<d<0 with the modified Cornacchia}
  { algorithm. Checks @x<>@y, but not p prime. Err<>0 if no solution exist.}
var
  a,b,t: mp_int;
  c2: integer;
label
  leave;
begin
  {Ref: [24] Cohen, Algorithm 1.5.3}
  {     [10] Crandall/Pomerance, Algorithm 2.3.13}

  {Error codes:
    -2:  d >= 0 or d <= -4p or p<2 or initial mp_error<>0
     4:  (p-x^2)/d is no square
     5:  (p-x^2)/d is no integer
   else  error code from mp_sqrtmod}

  Err := -2;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(d) or mp_not_init(p) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cornacchia4');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @x=@y then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_cornacchia4: @x=@y');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  mp_init3(a,b,t);
  if mp_error<>MP_OKAY then exit;

  {check -4p < d < 0, implemented as d<0 and |d|<|4p|}
  {set t=4p, will be used here and later for sqrt(4p)}
  mp_shl(p,2,t);
  c2 := mp_cmp_d(p,2);

  if (d.sign=MP_ZPOS) or (mp_cmp_mag(d,t)<>MP_LT) or (c2=MP_LT) then goto leave;

  if c2=MP_EQ then begin
    {handle case p=2}
    mp_add_d(d,8,t);
    if mp_is_square2(t,@x) then mp_set1(y);
    goto leave;
  end;

  {Step 1a: Solve b^2 = d mod p, Jacobi check (d|p) is done in mp_sqrtmod}
  mp_sqrtmod(d,p,b,Err);
  if Err=0 then begin
    {Step 1b: if b <> d mod 2 then replace b by p-b}
    if (b.pdigits^[0] and 1)<>(d.pdigits^[0] and 1) then mp_sub(p,b,b);
    {t := sqrt(4p)}
    mp_sqrt(t,t);
    {Euclid starting value a = 2p}
    mp_mul_2(p,a);
    {Step 2: Apply Euclidean algorithm to (a,b) until b <= sqrt(p)}
    while (mp_error=MP_OKAY) and mp_is_gt(b,t) do begin
      mp_mod(a,b,a);
      mp_exch(a,b);
    end;
    {Step 3: Set x=b, y = sqrt((4p-x^2)/|d|). If y is in N, then the}
    {solution is (x,y). Otherwise there is no solution, return Err}
    mp_sqr(b,t);
    mp_shl(p,2,a);
    mp_sub(a,t,a);
    mp_abs(d,t);
    mp_divrem(a,t,@a,@t);
    if mp_is0(t) then begin
      if mp_is_square2(a,@y) then
        mp_exch(x,b)
      else Err := 4;
    end
    else Err := 5;
  end;

leave:
  mp_clear3(a,b,t);
end;


{---------------------------------------------------------------------------}
procedure s_mp_cornacchia_ex(const a,b,p,q: mp_int; k: longint; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = m, m=p^k if k>0 otherwise m=p*q, a,b > 0, a+b < m,}
  { (a,m)=1 with Cornacchia's algorithm. Check @x<>@y, but not p,q prime.   }
  { Err <> 0 if error or no solution exist; normalisation x >= y if a=b.    }
var
  m,r,r2,s,t: mp_int;
  n: integer;
  aeqb: boolean;
label
  leave;
begin
  {Ref: A. Nitaj, L'algorithme de Cornacchia, Expositiones Mathematicae 13 (1995),}
  {     p.358-365; available from http://www.math.unicaen.fr/~nitaj/cornacchia.ps }
  {C.f. [24] Cohen, Algorithm 1.5.2,  [10] Crandall/Pomerance, Algorithm 2.3.12}

  {Error codes:
    -2:  a<1, b<1, a+b>m, k<0, or initial mp_error<>0
     5:  inverse of a mod p does not exists
     6:  (m-a*x^2)/b is no square
     7:  (m-a*x^2)/b is no integer
   else  error code from mp_sqrtmodpk/pq}

  Err := -2;
  if (mp_error<>MP_OKAY) or (k<0) then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(p) or mp_not_init(q) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_cornacchia_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @x=@y then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mp_cornacchia_ex: @x=@y');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {Test a<1 or b<1}
  if s_mp_is_le0(a) or s_mp_is_le0(b) then exit;

  mp_init5(m,r,r2,s,t);

  {Use p^2 if p and q are the same variable}
  if (k=0) and (@p=@q) then k:=2;

  if k>0 then mp_expt_int(p,k,m)
  else mp_mul(p,q,m);

  mp_add(a,b,t);
  if mp_is_gt(t,m) then goto leave;

  {here a>0, b>0, m>1. Remember if a=b for normalisation}
  aeqb := mp_is_eq(a,b);

  {Step 1a: Solve r^2 = -b/a mod m}
  if not mp_invmodf(a,m,t) then begin
    Err := 5;
    goto leave;
  end;
  mp_mulmod(t,b,m,t);
  s_mp_chs(t);
  mp_mod(t,m,t);
  if k>0 then mp_sqrtmodpk(t,p,k,r,Err)
  else mp_sqrtmodpq(t,p,q,r,r2,Err);

  n := 0;
  if Err=0 then begin
    repeat
      if n>0 then begin
        {Try second root if solving mod p*q and first root gives no solution}
        Err := 0;
        mp_exch(r,r2);
      end;
      {Step 1b: if r <= m/2, or equivalent 2r <= m, then replace r by m-r}
      mp_mul_2(r,t);
      if mp_is_le(t,m) then mp_sub(m,r,r);
      mp_div(m,a,t);
      mp_sqrt(t,t);
      mp_copy(m,s);
      {Step 2: Apply Euclidean algorithm to (r,s=m) until r <= t=sqrt(m/a)}
      while (mp_error=MP_OKAY) and mp_is_gt(r,t) do begin
        mp_mod(s,r,s);
        mp_exch(s,r);
      end;
      {Step 3: Set x=r, y = sqrt((m - a*r^2)/b). If y is in N, then the}
      {solution is (x,y). Otherwise there is no solution, return Err}
      mp_sqr(r,t);
      mp_mul(t,a,t);
      mp_sub(m,t,s);
      mp_divrem(s,b,@s,@t);
      if mp_is0(t) then begin
        if not mp_is_square2(s,@y) then Err := 6;
      end
      else Err := 7;
      inc(n);
    until (Err=0) or (k>0) or (n>1);
  end;
  if Err=0 then begin
    mp_exch(x,r);
    {Normalise solution: x >= y if a=b}
    if aeqb and mp_is_lt(x,y) then mp_exch(x,y);
  end;

leave:
  mp_clear5(m,r,r2,s,t);
end;


{---------------------------------------------------------------------------}
procedure mp_cornacchia_pk(const a,b,p: mp_int; k: longint; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = p^k, p prime; k,a,b > 0, a + b < p^k, (a,p)=1 with Cornacchia's  }
  { algorithm. Check @x<>@y, but not p prime. Err <> 0 if no solution exist. x >= y if a=b.}
begin
  if k<=0 then Err := -2
  else s_mp_cornacchia_ex(a,b,p,p,k,x,y,Err);
end;


{---------------------------------------------------------------------------}
procedure mp_cornacchia_pq(const a,b,p,q: mp_int; var x,y: mp_int; var Err: integer);
  {-Solve a*x^2 + b*y^2 = p*q, p,q prime; a,b > 0, a + b < p*q, (a,p*q)=1 with Cornacchia's  }
  { algorithm. Check @x<>@y, but not p,q prime. Err <> 0 if no solution exist. x >= y if a=b.}
begin
  s_mp_cornacchia_ex(a,b,p,q,0,x,y,Err);
end;


{---------------------------------------------------------------------------}
procedure mp_crt_setup(n: integer; const m: array of mp_int; var c: array of mp_int);
  {-Calculate CRT coefficients c[i] for pairwise co-prime moduli m[i], i=0..n-1.}
begin
  if not mp_crt_setupf(n,m,c) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_crt_setup: setup failed');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
  end;
end;


{---------------------------------------------------------------------------}
function mp_crt_setupf(n: integer; const m: array of mp_int; var c: array of mp_int): boolean;
  {-Calculate CRT coefficients c[i] for pairwise co-prime moduli m[i], i=0..n-1.}
  { Return true if c[i] are successfully calculated.}
label
  _err;
var
  u: mp_int;
  i,j: integer;
begin
  mp_crt_setupf := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init_multi(m) or mp_not_init_multi(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_crt_setupf');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  dec(n);
  if (n > high(m)) or (n > high(c)) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_crt_setupf: invalid n');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_init(u);
  if mp_error<>MP_OKAY then exit;
  {Step 1 of Garner's algorithm for CRT, HAC Alg. 14.71}
  for i:=0 to n do begin
    mp_set1(c[i]);
    for j:=0 to pred(i) do begin
      if (mp_error<>MP_OKAY) or not mp_invmodf(m[j],m[i],u) then goto _err;
      mp_mulmod(u,c[i],m[i],c[i]);
    end;
  end;
  mp_crt_setupf := mp_error=MP_OKAY;

_err:
  mp_clear(u);
end;


{---------------------------------------------------------------------------}
procedure mp_crt_single(n: integer; const m,v: array of mp_int; var x: mp_int);
  {-Calculate x with x mod m[i] = v[i], pairwise co-prime moduli m[i], i=0..n-1.}
  { Use mp_crt_setup/calc if more than one system shall be solved}
var
  u,p,ci: mp_int;
  i: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init_multi(m) or mp_not_init_multi(v) or mp_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_crt_single');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  dec(n);
  if (n<0) or (n>high(m)) or (n>high(v)) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_crt_single: invalid n');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {start with x = v[0] mod m[0]}
  mp_mod(v[0],m[0],x);
  if n=0 then exit; {trivial case, don't use local mp_ints}

  mp_init3(u,p,ci);
  if mp_error=MP_OKAY then begin
    {already computed: mp_mod(v[0],m[0],x)}
    mp_set1(p);
    for i:=1 to n do begin
      {p = product(m[j], j=0..i-1)}
      mp_mul(p,m[i-1],p);
      {calculate c[i] from mp_crt_setup 'on the fly'}
      mp_invmod(p,m[i],ci);
      mp_submod(v[i],x,m[i],u);
      mp_mulmod(u,ci,m[i],u);
      mp_mul(u,p,u);
      mp_add(x,u,x);
    end;
    mp_clear3(u,p,ci);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_crt_solve(n: integer; const m,c,v: array of mp_int; var x: mp_int);
  {-Calculate x with x mod m[i] = v[i], pairwise co-prime moduli m[i], i=0..n-1.}
  { Coefficients c[i] must be precalculated with mp_crt_setup}
var
  u,p: mp_int;
  i: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init_multi(m) or mp_not_init_multi(c) or mp_not_init_multi(v) or mp_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_crt_solve');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  dec(n);
  if (n<0) or (n>high(m)) or (n>high(c)) or (n>high(v)) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_crt_solve: invalid n');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {start with x = v[0] mod m[0]}
  mp_mod(v[0],m[0],x);
  if n=0 then exit; {trivial case, don't use local mp_ints}

  mp_init2(u,p);
  if mp_error=MP_OKAY then begin
    {Step 2 of Garner's algorithm for CRT, HAC Alg. 14.71}
    {already computed: mp_mod(v[0],m[0],x)}
    {p accumulates the product of the m[i]}
    mp_set1(p);
    {Step 3}
    for i:=1 to n do begin
      {p = product(m[j], j=0..i-1)}
      mp_mul(p,m[i-1],p);
      mp_submod(v[i],x,m[i],u);
      mp_mulmod(u,c[i],m[i],u);
      mp_mul(u,p,u);
      mp_add(x,u,x);
    end;
    mp_clear2(u,p);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_dfact(n: longint; var a: mp_int);
  {-Calculate double factorial a = n!!, a=n*(n-2)..., lowest term 1 or 2; a=0 for n < -3}
  { Note: for protected mode code n=$FFFF is near the maximum value, but}
  { for real mode BP7 one might get a heap overflow even for smaller n!}
const
  {double factorial for small arguments}
  sdfact: array[-3..19] of longint =
            (-1,0,1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,
             645120,2027025,10321920,34459425,185794560,654729075);
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_dfact');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n<20 then begin
    if n<-3 then mp_zero(a)
    else mp_set_int(a,sdfact[n])
  end
  else if odd(n) then mp_OddProd(0,n shr 1,a)
  else begin
    { (2n)!! = 2^n * n!}
    n := n shr 1;
    mp_fact(n,a);
    mp_shl(a,n,a);
  end;
end;


{---------------------------------------------------------------------------}
const
  {factorial for small arguments}
  sm_fact: array[0..12] of longint = (1,1,2,6,24,120,720,5040,40320,362880,
                                      3628800,39916800,479001600);

{---------------------------------------------------------------------------}
procedure s_mp_recsplit_fact(n: word; var a: mp_int);
  {-Calculate a = factorial(n) using Recursive Split, error if n > MaxFact}
const
  MAXL = 14;
var
  np : word;                     {Luschny's "private long N"}
  st : array[0..MAXL] of mp_int; {local 'stack', avoids init/clear in product}
  lev: integer;                  {stack level >=2; 0,1 for local mp_ints}

  procedure product(var p: mp_int; l: word);
    {-Recursive product of odd numbers (depends on np and calling sequence)}
  const
    lmin={$ifdef BIT16}8{$else}32{$endif}; {min value of l for recursion}
  begin
    if mp_error<>MP_OKAY then exit;
    if l<lmin then begin
      inc(np,2);
      mp_set_w(p, np);
      {use while, this allows variable limits}
      while l>1 do begin
        inc(np, 2);
        mp_mul_w(p, np, p);
        dec(l);
      end;
    end
    else begin
      if lev>=MAXL then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXRange.Create('s_mp_recsplit_fact: lev out of range');
          {$else}
            RunError(MP_RTE_RANGE);
          {$endif}
        {$else}
          set_mp_error(MP_RANGE);
          exit;
        {$endif}
      end;
      inc(lev);
      product(p, l shr 1);
      if mp_error=MP_OKAY then begin
        product(st[lev], l-(l shr 1));
        mp_mul(p, st[lev], p);
      end;
      dec(lev);
    end;
  end;

var
  i,high,len: word;

begin
  {Based of P.Luschny's FactorialSplitRecursive.java, available at}
  {http://www.luschny.de/math/factorial/FastFactorialFunctions.htm}
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_recsplit_fact');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n <= 12 then mp_set_int(a, sm_fact[n])
  else begin
    {initialize "private long N"}
    np := 1;
    {Create and initialize local mp_ints. Note that using the stack is an}
    {improvement only for WIN32 (up to 30%) where memory allocation seems}
    {to be a bit time consuming. No improvement for BP7 and FPC/go32v2.}
    mp_init_multi(st);
    if mp_error<>MP_OKAY then exit;
    mp_set(st[0], 1);
    mp_set(st[1], 1);

    high := 1;
    for i:=bitsize32(n)-2 downto 0 do begin
      len  := high;
      high := n shr i;
      if not odd(high) then dec(high);
      len := (high-len) shr 1;
      if len>0 then begin
        {initialize stack level, product uses levels lev+1,...}
        lev := 1;
        {use a for product result}
        product(a,len);
        mp_mul(st[0],a,st[0]);
        mp_mul(st[1],st[0],st[1]);
      end;
    end;
    {Luschny's final shift is n - popcount(n)}
    mp_shl(st[1], n - popcount16(n), a);
    mp_clear_multi(st);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_borsch_fact(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) with Borwein/Schoenhage prime factorization method}
{$ifdef BIT16}
const
  FMAXPRIME = $F000 div 4;
{$else}
const
  FMAXPRIME = 8000000;
{$endif}
type
  TFTab = array[0..FMAXPrime] of longint;
var
  primtab,factab,exptab: ^TFTab;
  expo,maxe,k,p,nh,i,bs,maxi,tsize: longint;
  prod: mp_int;
  sieve: TSieve;
begin

  {Compute the factorial n! based on the prime factorization of n. Some of the }
  {basic facts are described in P.B. Borwein, On the complexity of calculating }
  {factorials, Journal of Algorithms, Vol.6, pp. 376-380, 1985. Available as   }
  {http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P29.pdf. See also P. Luschny}
  {http://www.luschny.de/math/factorial/java/FactorialPrimeSchoenhage.java.html}
  {http://www.luschny.de/math/factorial/java/FactorialPrimeBorwein.java.html   }
  {or http://en.literateprograms.org/Factorials_with_prime_factorization_(C)   }

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_borsch_fact');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n <= 12 then begin
    mp_set_int(a, sm_fact[n]);
    exit;
  end;

  nh    := n div 2;
  maxe  := 1;
  maxi  := primepi32(n);
  tsize := (maxi+8)*sizeof(longint);

  primtab := mp_getmem(tsize);
  mp_set1(a);

  {Step 1 of Borwein's method: construct primes}
  if n <= $FFFF then begin
    {Use static 16-bit primes}
    for i:=1 to maxi do primtab^[i] := primes16[i];
  end
  else begin
    {Use prime generator if n has more than 16 bits}
    prime_sieve_init(sieve, 2);
    for i:=1 to maxi do primtab^[i] := prime_sieve_next(sieve);
    prime_sieve_clear(sieve);
  end;

  {Step 2: Building the exponents table for the primes <= n:  We skip}
  {p=2 and start with p=3, the powers of 2 are shifted in at the end.}
  {maxe is the maximum exponent of the primes, used to avoid loops   }
  {and useless squarings in the second (accumulation) part.          }

  factab  := mp_getmem(tsize);
  exptab  := mp_getmem(tsize);

  for i:=2 to maxi do begin
    p := primtab^[i];
    if p>nh then expo := 1
    else begin
      k := n div p;
      expo := k;
      while k>0 do begin
        k := k div p;
        inc(expo,k);
      end;
    end;
    exptab^[i] := expo;
    if expo>maxe then maxe := expo;
  end;

  mp_init(prod);
  bs := bitsize32(n);
  while maxe shr bs = 0 do dec(bs);

  {Combined steps 3 and 4: Compute and accumulate the alphas = prod}
  while (bs >= 0) and (mp_error = MP_OKAY) do begin
    k := 0;
    for i:=2 to maxi do begin
      if odd(exptab^[i] shr bs) then begin
        factab^[k] := primtab^[i];
        inc(k);
      end;
    end;
    mp_sqr(a,a);
    if k>0 then begin
      mp_prod_int(prod,factab^,k);
      mp_mul(a,prod,a);
    end;
    dec(bs);
  end;

  {Multiply by powers of two}
  mp_shl(a,n - popcount32(n),a);

  mp_clear(prod);
  mp_freemem(pointer(primtab), tsize);
  mp_freemem(pointer(factab),  tsize);
  mp_freemem(pointer(exptab),  tsize);
end;


{---------------------------------------------------------------------------}
procedure mp_fact(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) using Recursive Split or Borwein/Schoenhage}
  { prime factorization method, error if n > MaxFact}
begin
  if (n < 0) or (n > MaxFact) then begin
    {$ifdef MPC_ArgCheck}
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_fact: invalid n');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
      {$endif}
    {$endif}
    exit;
  end
  else if n<1000 then s_mp_recsplit_fact(n,a)
  else s_mp_borsch_fact(n,a);
end;


{---------------------------------------------------------------------------}
procedure mp_fermat(n: word; var fn: mp_int);
  {-Return nth Fermat number, fn = 2^(2^n)+1 (MP_RANGE error for n>MaxFermat)}
begin
  if mp_error<>MP_OKAY then exit;
  {Check n because 1 shl n is calculated as 1 shl (n mod 32) for n>31}
  if n>MaxFermat  then begin
    {$ifdef MPC_HaltOnArgCheck}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_fermat: n too large');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {Init check in mp_2expt}
  mp_2expt(fn, longint(1) shl n);
  if mp_error=MP_OKAY then begin
    fn.pdigits^[0] := fn.pdigits^[0] or 1;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_fib2(n: longint; var fn,f1: mp_int);
  {-Calculate two Fibonacci numbers fn=fib(n), f1=fib(n-1), n>=0}
var
  t: mp_int;
  bk: longint;
{$ifdef MPC_ArgCheck}
const
  toobig = trunc(MaxMersenne*1.5);
{$endif}
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(fn) or mp_not_init(f1) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_fib2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if n<0 then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_fib2: n<0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
    if n>toobig then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_fib2: n too large');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
    if @fn=@f1 then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_fib2: @fn=@f1');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  {initialize f1=0, fn=1}
  mp_zero(f1);
  mp_set1(fn);
  {easy out for n=0,1}
  if n=0 then mp_exch(fn,f1);
  if n<2 then exit;

  {create temporary variable}
  mp_init(t);

  {fn=F[n] and f1=F[n-1] are calculated simultaneously using a}
  {"left-to-right" binary algorithm from high to low bits of n}
  {   F[k+1] =   F[k]   + F[k-1],  F[0]=1, F[1]=1             }
  {  F[2k+1] = 4*F[k]^2 - F[k-1]^2 + 2*(-1)^k                 }
  {  F[2k-1] =   F[k]^2 + F[k-1]^2                            }
  {  F[2k]   = F[2k+1]  - F[2k-1]                             }

  {get highest bit of n, loop terminates because n>0}
  bk := longint($40000000);
  while n and bk = 0 do bk := bk shr 1;

  while bk>1 do begin
    {bk is the highest unprocessed bit}
    {f1 = F[k-1]^2}
    mp_sqr(f1,f1);
    {t  = F[k]^2  }
    mp_sqr(fn,t);

    {F[2k+1] = 4*F[k]^2 - F[k-1]^2 + 2*(-1)^k}
    {     fn =    4*t   -     f1   + 2*(-1)^k}
    mp_shl(t,2,fn);
    mp_sub(fn,f1,fn);
    {add 2*(-1)^k, note: mp_digits are positive so use add_d or sub_d}
    if n and bk = 0 then mp_add_d(fn,2,fn) else mp_sub_d(fn,2,fn);

    {F[2k-1] = F[k]^2 + F[k-1]^2}
    {     f1 =    t   +     f1  }
    mp_add(t,f1,f1);

    {Get next lower bit of n. If it is 0 then F[2k] and F[2k-1]}
    {are used, otherwise F[2k+1] and F[2k] are used}
    bk := bk shr 1;
    if n and bk = 0 then begin
      {f1 = F[2k-1], fn = F[2k]  }
      mp_sub(fn,f1,fn);
    end
    else begin
      {f1 = F[2k]),  fn = F[2k+1]}
      mp_sub(fn,f1,f1);
    end;
  end;
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_fib(n: longint; var fn: mp_int);
  {-Calculate Fibonacci number fn=fib(n), fib(-n)=(-1)^(n-1)*fib(n)}
var
  f1: mp_int;
const
  smfib: array[0..24] of word = (0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,
                                 1597,2584,4181,6765,10946,17711,28657,46368);
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(fn) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_fib');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if abs(n)<25 then mp_set_w(fn,smfib[abs(n)])
  else begin
    {create temporary variable for fib(n-1)}
    mp_init(f1);
    if mp_error=MP_OKAY then begin
      mp_fib2(abs(n),fn,f1);
      mp_clear(f1);
    end;
  end;
  {adjust sign if n negative}
  if (n<0) and not odd(n) then mp_chs(fn,fn);
end;


{---------------------------------------------------------------------------}
function mp_isMersennePrime(p: longint): boolean;
  {-Lucas-Lehmer test for Mersenne number m=2^p-1, HAC 4.37}
var
  k,n,kmax,t: longint;
  r8 : integer;
  m,u: mp_int;
begin
  mp_isMersennePrime := false;

  {If p not prime, then m not prime}
  if (p<2) or not IsPrime32(p) then exit;

  {Here p is prime, check for p=2,3 or p a Sophie Germain prime}
  {Lucas-Lehmer test assumes p>=3}
  if p<=3 then begin
    mp_isMersennePrime := true;
    exit;
  end;

  {Ribenboim, Ch.2, VII: if p=3 (mod 4), p>3, q=2p+1 prime then q is a}
  {factor of m. This skips about 6% - 10% full tests for p < 100000000}
  if (p and 3 = 3) and (p<MaxLongint) and IsPrime32(succ(p shl 1)) then exit;

  mp_init2(m,u);
  if mp_error<>MP_OKAY then exit;

  {calculate 2^p-1}
  mp_mersenne(p,m);

  {Factors n of m have form n=1,7 mod 8 and n=1 mod 2p, ie n=2kp+1}
  {see http://www.utm.edu/research/primes/notes/proofs/MerDiv.html}
  {Search factor 2kp+1 with 1 <= k < min(MaxLongint/p, 2*p+1)}

  {smallest possible factor}
  n := succ(2*p);
  kmax := MaxLongint div p;
  if n<kmax then kmax := n;
  if p<60 then begin
    t := (1 shl ((p+1) div 2)) div p;
    if kmax>t then kmax := t;
  end;

  for k:=1 to kmax do begin
    r8 := n and 7;
    if ((r8=1) or (r8=7)) and IsPrime32(n) then begin
      mp_mod_int(m,n,t);
      if t=0 then begin
        {Check for spurious factor for small p}
        if mp_is_longint(m,t) and (n=t) then break;
        mp_clear2(m,u);
        exit;
      end;
    end;
    inc(n,2*p);
  end;

  {run Lucas-Lehmer test}
  mp_set(u,4);
  {Use unrestricted diminished radix reduction. Here we know that}
  {m>0 and the setup value is d=1, so skip test/setup functions  }
  for k:=1 to p-2 do begin
    {calculate u := u^2 - 2 mod m}
    mp_sqr(u,u);
    mp_sub_d(u,2,u);
    {if u<0 make u positive, no reduction}
    if u.sign=MP_NEG then mp_add(u,m,u)
    else mp_reduce_2k(u,m,1);
    if mp_error<>MP_OKAY then break;
  end;
  mp_isMersennePrime := mp_iszero(u) and (mp_error=MP_OKAY);
  mp_clear2(m,u);
end;


{---------------------------------------------------------------------------}
function s_mp_is_pth_power(const a: mp_int; p: longint; var r: mp_int): boolean;
  {-Return true if a is pth power, then a=r^p. a>0, p>2 prime, no checks}
var
  amod64: mp_digit;
  w: word;
  j: integer;
  c,q,t: longint;
const
  {const arrays calculated with t_pmtabs.pas}
  ba_03_13: array[0..1] of byte = ($23,$11);
  ba_03_61: array[0..7] of byte = ($0b,$0b,$90,$19,$66,$02,$34,$14);
  ba_03_63: array[0..7] of byte = ($03,$01,$00,$18,$18,$00,$80,$40);
  ba_03_64: array[0..7] of byte = ($ab,$ab,$aa,$ab,$aa,$ab,$aa,$ab);

  ba_05_25: array[0..3] of byte = ($83,$00,$04,$01);
  ba_05_41: array[0..5] of byte = ($0b,$42,$00,$08,$41,$01);
  ba_05_44: array[0..5] of byte = ($03,$18,$a0,$00,$03,$08);
  ba_05_64: array[0..7] of byte = ($ab,$aa,$aa,$aa,$ab,$aa,$aa,$aa);

  ba_07_29: array[0..3] of byte = ($03,$10,$02,$10);
  ba_07_43: array[0..5] of byte = ($c3,$00,$00,$00,$30,$04);
  ba_07_49: array[0..6] of byte = ($03,$00,$0c,$c0,$00,$00,$01);

  ba_11_23: array[0..2] of byte = ($03,$00,$40);
  ba_11_36: array[0..4] of byte = ($b3,$2b,$9b,$ba,$09);

  ba_13_53: array[0..6] of byte = ($03,$00,$80,$40,$00,$00,$10);

  ba_xx_64: array[0..7] of byte = ($ab,$aa,$aa,$aa,$aa,$aa,$aa,$aa);

  {----------------------------------------------------------}
  function isclr(const ba: array of byte; k: word): boolean;
    {-Test is bit k is zero}
  var
    i: integer;
  begin
    i := k shr 3;
    isclr := (i > high(ba)) or (ba[i] and (1 shl (k and 7)) = 0);
  end;

begin
  s_mp_is_pth_power := false;

  if (mp_error=MP_OKAY) and (a.used>0) and isprime32(p) then begin

    amod64 := a.pdigits^[0] and 63;

    if p=3 then begin
      if isclr(ba_03_64, amod64) then exit;
      s_mp_mod_w(a,13*61*63,w);
      if isclr(ba_03_13, w mod 13) or isclr(ba_03_61, w mod 61) or isclr(ba_03_63, w mod 63) then exit;
    end
    else if p=5 then begin
      if isclr(ba_05_64, amod64) then exit;
      s_mp_mod_w(a,25*41*44,w);
      if isclr(ba_05_25, w mod 25) or isclr(ba_05_41, w mod 41) or isclr(ba_05_44, w mod 44) then exit;
    end
    else begin
      if isclr(ba_xx_64, amod64) then exit;
      if p=7 then begin
        s_mp_mod_w(a,29*43*49,w);
        if isclr(ba_07_29, w mod 29) or isclr(ba_07_43, w mod 43) or isclr(ba_07_49, w mod 49) then exit;
      end
      else if p=11 then begin
        s_mp_mod_w(a,23*36,w);
        if isclr(ba_11_23, w mod 23) or isclr(ba_11_36, w mod 36) then exit;
      end
      else if p=13 then begin
        s_mp_mod_w(a,53,w);
        if isclr(ba_13_53, w) then exit;
      end;
    end;

    {Check if a is no pth power residue of m, (m=p^2 or m=q=t*p+1, q prime) }
    {Theorem 6.18 from A.Adler, J.E.Coury [26]: If m>0 has a primitive root,}
    {gcd(a,p)=1, then x^p = a (mod m) has a solution if and only if         }
    {a^(phi(m)/gcd(p,phi(m))) = 1 mod m.  Since m=p^2 or m=q are suitable,  }
    {we know that if  a<>0 mod m  and  a^(phi(m)/gcd(p,phi(m))) <> 1 mod m  }
    {then x^p = a has no solutions mod m, and x^p = a has no solution in Z. }
    {Since p is prime, gcd(a,p)=1 if a mod p <> 0}

    {First tests only if p^2 < MaxLongint}
    if p<=46340 then begin
      c := sqr(p);
      s_mp_mod_w(a,word(p),w);
      mp_mod_int(a,c,t);
      {w = a mod p; t = a mod p^2}
      if w<>0 then begin
        {a is pth power residue mod p^2 if a mod p <> 0 and}
        {1 = a^(phi(p^2)/gcd(p,phi(p^2))) = a^(p-1) mod p^2}
        if exptmod32(t,pred(p),c)<>1 then exit;
      end
      else begin
        {Here a mod p = 0, i.e. a = x*p. Since p>2, a must be a multiple}
        {of p^3 if it is a pth power. If  a mod p^2 <> 0, a is no power!}
        if t<>0 then exit
        else if p<=1290 then begin
          {Check a mod p^3 = 0}
          mp_mod_int(a,c*p,t);
          if t<>0 then exit;
        end;
      end;
    end;

    {get primes q = t*p + 1}
    q := 1;
    t := 0;
    {Note: Loop count for 'sieve tests' approximately equal to the value }
    {2*log2(log2(a)) / log2(p). Ref: E.Bach, J.Sorenson: Sieve Algorithms}
    {for Perfect Power Testing, Algorithmica 9 (1993), 313-328}
    for j:=0 to (2*bitsize32(mp_bitsize(a))) div bitsize32(p) do begin
      {find next prime q = t*p + 1}
      c := p+p;
      repeat
        inc(q,c);
        inc(t,2);
      until (q<0) or isprime32(q);
      {break if 32-bit primes exceeded}
      if q<2 then break;
      mp_mod_int(a,q,c);
      {Check a^(phi(q)/gcd(p,phi(q))) <> 1 mod q}
      {phi(q)=q-1, gcd(p,phi(q)) = gcd(p,p*t) = p}
      {phi(q)/gcd(p,phi(q)) = p*t/p = t}
      if c<>0 then begin
        if exptmod32(c,t,q)<>1 then begin
          {a is no pth power residue mod q}
          exit;
        end;
      end
      else begin
        {a mod q = 0, no power if a mod q^2 <> 0}
        if q<=46340 then begin
          mp_mod_int(a,q*q,c);
          if c<>0 then exit;
        end;
      end;
    end;
  end;
  {a has survived all fast checks, compute r=floor(a^(1/p)) and check r^p=a}
  s_mp_is_pth_power := s_mp_n_root2(a,p,r,nil) and (mp_error=MP_OKAY);
end;


{---------------------------------------------------------------------------}
function mp_is_pth_power(const a: mp_int; p: longint; var r: mp_int): boolean;
  {-Return true if a is pth power, a>0, p prime. If true, calculate r with a=r^p}
begin
  if p=2 then begin
    mp_is_pth_power := mp_is_square2(a,@r);
    exit;
  end;
  mp_is_pth_power := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(r) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_pth_power');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if a.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_is_pth_power: a<0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  if (p<2) or not isprime32(p) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_is_pth_power: p not prime');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  mp_is_pth_power := s_mp_is_pth_power(a,p,r);
end;


{---------------------------------------------------------------------------}
procedure mp_is_power(const a: mp_int; var b: mp_int; var p: longint);
  {-Calculate smallest prime p with a=b^p; p=1,b=a if a is no power}
var
  r,t,s: mp_int;
  i,l,n: longint;
  k: integer;
  f,fmax: mp_digit;
  sa: word;
  amod8: mp_digit;
  FR: TFactors32;
  fac: double;

{Constants for trial division part, c.f. the technical report   }
{E.Bach, J.Sorensen, Sieve Algorithms for Perfect Power Testing,}
{from http://research.cs.wisc.edu/techreports/1989/TR852.pdf    }

{p=pm[i] is roughly the smallest prime with p*log2(p)^2 >= log2(a), index i}
{is bitsize(bitsize(a)), pf[i] is slightly >= ln(2)/ln(nextprime(pm[i]+1))}
const
  pm: array[9..21] of mp_digit=
       (23, 41, 61, 97, 157, 257, 431, 727, 1249, 2143, 3727, 6547, 11527);
  pf: array[9..21] of single =
       (0.205847, 0.184289, 0.164851, 0.150191, 0.136078, 0.124395, 0.114179,
        0.105068, 0.0971057, 0.0903169, 0.0842736, 0.0788800, 0.0740989);

  {--------------------------------------------------}
  procedure processprime(d: mp_digit);
    {-d is a prime factor of t: t = d^n*s, process s with prime factorization of n}
    { if p<>1 then |a| = r^p}
  var
    q: longint;
    j: integer;
  begin
    mp_valrem(t,d,s,n);
    {t = d^n*s is a power if s=1 or s=r^p with p|n}
    PrimeFactor32(n,FR);
    with FR do begin
      for j:=1 to pcount do begin
        q := primes[j];
        if (q=2) and (sa=MP_NEG) then continue;
        if mp_is1(s) then begin
          p := q;
          mp_set_pow(r,d,n div p);
          exit;
        end
        else if s_mp_is_pth_power(s,q,r) then begin
          p := q;
          mp_set_pow(s,d,n div p);
          mp_mul(s,r,r);
          exit;
        end;
      end;
    end;
  end;

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_power');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {remember sign(a) and |a| mod 8}
  sa := a.SIGN;
  amod8 := a.pdigits^[0] and 7;

  {init default answer}
  p := 1;

  {check trivial cases}
  l := mp_bitsize(a);
  if (amod8 and 3 = 2) or (l<4) then begin
    {b=a^1 for a=-7..7; except 4=2^2. Further b=a^1 for a=2 mod 4}
    {because a^p <> 2 mod 4 for all primes}
    if (amod8=4) and (sa=MP_ZPOS) then begin
      {a=4, return b=2, p=2}
      mp_set(b,2);
      p := 2;
    end
    else mp_copy(a,b);
    exit;
  end;

  {if a is positive, check if a is square}
  if (sa=MP_ZPOS) and mp_is_square2(a,@b) then begin
    p := 2;
    exit;
  end;

  mp_init3(r,s,t);
  {perform calculation for t=abs(a)}
  mp_abs(a,t);

  if mp_isodd(t) then begin
    {fac = ln(2)/ln(p1) with p1 first prime >= 3 for worst case}
    if l<30 then begin
      f := 0;
      fac := 0.63093;
    end
    else begin
      k := bitsize32(l);
      if k>21 then k := 21;
      if k<9  then k := 9;
      {$ifdef MP_16BIT}
        while pm[k]>MP_DIGIT_MAX do dec(k);
      {$endif}
      fmax := pm[k];
      fac  := pf[k];
      mp_small_factor(t,3,fmax,f);
    end;
    if f<>0 then processprime(f)
    else begin
      {No small factor <= fmax found, set l = maximum prime to test}
      {for the worst case t = nextprime(fmax+1)^p. Because t is no }
      {square the loop through the i'th prime powers starts with 3.}
      l := 1 + trunc(fac*l);
      i := 3;
      while (i<=l) and (mp_error=MP_OKAY) do begin
        if s_mp_is_pth_power(t,i,r) then begin
          if not mp_is1(r) then p := i;
          break;
        end;
        i := nextprime32(i+1);
      end;
    end;
  end
  else begin
    {t is even}
    if (mp_popcount(t)=1) and isprime32(l-1) then begin
      {t=2^p with p=l-1}
      p := l-1;
      mp_set(r,2);
    end
    else begin
      {a = 2^n * t, perform optimized processprime(2)}
      mp_makeodd(t,t,n);
      PrimeFactor32(n,FR);
      with FR do begin
        for k:=1 to pcount do begin
          i := primes[k];
          if (i=2) and (sa=MP_NEG) then continue;
          if s_mp_is_pth_power(t,i,r) then begin
            p := i;
            mp_shl(r,n div p, r);
            break;
          end;
        end;
      end;
    end;
  end;

  {if no p found set b=a}
  if p=1 then begin
    {if no p>1 found set b=a}
    mp_copy(a,b)
  end
  else if p=2 then begin
    {Bugfix 0.9.01: if p=2 and a<0 then set p=1 and b=a. The only case}
    {should be a=-4 resulting from |a|=2^p with p even prime!}
    if sa=MP_NEG then begin
      p := 1;
      mp_copy(a,b)
    end
    else mp_exch(r,b);
  end
  else begin
    {set b=sign(a)*r, (note p is odd)}
    if sa=MP_NEG then mp_chs(r,b)
    else mp_exch(r,b);
  end;
  mp_clear3(r,s,t);
end;


{---------------------------------------------------------------------------}
procedure mp_is_power_max(const a: mp_int; var b: mp_int; var k: longint);
  {-Calculate largest k with a=b^k; k=1,b=a if a is no power}
var
  n: longint;
begin
  {Arg checking is done in mp_copy}
  mp_copy(a,b);
  if mp_error<>MP_OKAY then exit;

  {invariant a=b^k}
  k := 1;
  {done if |a|=1}
  if mp_is1a(b) then exit;

  n := 1;
  repeat
    {set b'=b^n, k'=k*n until n=1}
    mp_is_power(b,b,n);
    k := k*n;
  until (n=1) or (mp_error<>MP_OKAY);
end;


{---------------------------------------------------------------------------}
procedure mp_nextpower_ex(const a: mp_int; var b,c: mp_int; var p: longint);
  {-Return smallest perfect power b := c^p >= a}
var
  m,d,x,t: mp_int;
  k,l: longint;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_nextpower');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Checked against https://oeis.org/a001597/b001597.txt}
  if mp_cmp_d(a,2)=MP_LT then begin
    mp_set(b,1);
    mp_set(c,1);
    p := 1;
    exit;
  end;

  if mp_is_longint(a,k) and (k<17) then begin
    if k<=4 then begin
      mp_set(b,4);
      mp_set(c,2);
      p := 2;
    end
    else if k<=8 then begin
      mp_set(b,8);
      mp_set(c,2);
      p := 3;
    end
    else if k=9 then begin
      mp_set(b,9);
      mp_set(c,3);
      p := 2;
    end
    else begin
      mp_set(b,16);
      mp_set(c,4);
      p := 2;
    end;
    exit;
  end;

{
  mp_is_power(a,c,p);
  if p>1 then begin
    mp_expt_int(c,p,b);
    exit;
  end;
}

  mp_init4(d,m,t,x);
  l := mp_bitsize(a)+1;
  mp_mul_2k(a,1,m);
  k := 2;
  while k <= l do begin
    mp_n_root2(a,k,t,@d);

    {Done if remainder or root is zero}
    if d.used=0 then begin
      mp_exch(t,x);
      p := k;
      break;
    end;
    if t.used=0 then break;

    mp_inc(t);
    mp_expt_int(t,k,d);
    mp_sub(d,a,d);
    if mp_is_lt(d,m) then begin
      mp_exch(m,d);
      p := k;
      mp_exch(t,x);
      if (m.used=0) or (mp_is1(m)) then break;
    end;
    k := nextprime32(k+1);
  end;
  mp_copy(x,c);
  mp_expt_int(x,p,b);
  mp_clear4(d,m,t,x);
end;


{---------------------------------------------------------------------------}
procedure mp_nextpower(var a: mp_int);
  {-Return smallest perfect power >= a}
var
  b,c: mp_int;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  mp_init2(b,c);
  if mp_error=MP_OKAY then begin
    mp_nextpower_ex(a,b,c,p);
    mp_exch(a,b);
    mp_clear2(b,c);
  end;
end;

{---------------------------------------------------------------------------}
function s_mp_is_pprime_ex(const a: mp_int; smin,smax: mp_digit): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32); trial division}
  { from smin to smax; no init check, no check if a is less than 2<31 }
var
  f: mp_digit;
  n: longint;
const
  DC_MAXLONG = (30+DIGIT_BIT) div DIGIT_BIT;
begin
  if (a.used<=DC_MAXLONG) and mp_is_longint(a,n) then begin
    s_mp_is_pprime_ex := IsPrime32(n);
    exit;
  end;
  {Default to not prime}
  s_mp_is_pprime_ex := false;

  {check if a has small factor up to min(MP_DIGIT_MAX,smax)}
  mp_small_factor(a, smin, smax, f);
  if (mp_error<>MP_OKAY) or (f<>0) then exit;

  {Test for BPSW (Baillie-Pomerance-Selfridge-Wagstaff) probable prime:}
  {First check if a is spsp(2), exit if not}
  if not mp_is_spsp_d(a,2) then exit;

  {then do a Strong Lucas probable prime test}
  s_mp_is_pprime_ex := s_mp_is_lpsp(a, true) and (mp_error=MP_OKAY);

end;


{---------------------------------------------------------------------------}
function mp_is_pprime_ex(const a: mp_int; smax: mp_digit): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32); trial division up to smax}
var
  n: longint;
begin
  {Default to not prime}
  mp_is_pprime_ex := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_pprime_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Negative numbers and 0 are not prime}
  if s_mp_is_le0(a) then exit;

  {uses IsPrime32 for less than 32 bit}
  if mp_is_longint(a,n) then begin
    mp_is_pprime_ex := (mp_error=MP_OKAY) and (n>0) and IsPrime32(n);
    exit;
  end;

  {here we have a with more than one mp_digit and at least 32 bit}
  mp_is_pprime_ex := s_mp_is_pprime_ex(a, 2, smax and MP_DIGIT_MAX);

end;


{---------------------------------------------------------------------------}
function mp_is_pprime(const a: mp_int): boolean;
  {-Test if a is prime (BPSW probable prime if a>2^32)}
  { Trial division up to min(MP_DIGIT_MAX,$3FFF)}
const
  smax = mp_digit(MP_DIGIT_MAX and $3FFF);  {max. for small factor}
begin
  mp_is_pprime := mp_is_pprime_ex(a, smax);
end;


{---------------------------------------------------------------------------}
function mp_is_primepower(const a: mp_int; var b: mp_int; var k: longint): boolean;
  {-Return true if a=b^k, b prime, k>1, otherwise false and a=b^k, k=1 if no power}
var
  n: longint;
begin
  mp_is_primepower := false;

  {invariant a=b^k}
  mp_copy(a,b);
  k := 1;

  {Arg checking is done in mp_copy}
  if mp_error<>MP_OKAY then exit;

  {Numbers <= 1 are no prime powers}
  if s_mp_is_le0(b) or mp_is1(b) then exit;

  {Handle case of even a}
  if mp_iseven(b) then begin
    if mp_is_pow2(b,k) then begin
      mp_set(b,2);
      mp_is_primepower := (mp_error=MP_OKAY) and (k>1);
    end;
  end
  else begin
    n := 1;
    repeat
      {set b'=b^n, k'=k*n until n=1}
      mp_is_power(b,b,n);
      k := k*n;
    until (n=1) or (mp_error<>MP_OKAY);
    {a=b^k with k>1. Check if b is prime}
    if (k>1) and mp_is_pprime(b) then begin
      mp_is_primepower := mp_error=MP_OKAY;
      exit;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_is_slpsp_alt(const a: mp_int): boolean;
  {-Strong Lucas probable prime test for a. Lucas test is }
  { done for the first p=2k+1 with mp_jacobi(p^2-4,a) = -1}
const
  bmax = $7FFF; {max p for Jabobi test)}
  bsqr = 127;   {p index for perfect square test, must be odd}

label
  leave;

var
  b: integer;
  d,i: longint;
  mu,p,t,v: mp_int;
begin

  {initialize function result with false}
  mp_is_slpsp_alt := false;
  if mp_error<>MP_OKAY then exit;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_slpsp_alt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {easy outs: for a<3 return true only if a=2}
  b := mp_cmp_d(a,2);
  if (a.sign=MP_NEG) or (b=MP_LT) then exit;
  if mp_iseven(a) then begin
    mp_is_slpsp_alt := (b=MP_EQ) and (mp_error=MP_OKAY);
    exit;
  end;

  {create temporary mp_ints}
  mp_init4(mu,p,t,v);
  if mp_error<>0 then exit;

  {generate sequence without squaring: p(k)=2k+1, D(k)=p(k)^2 - 4}
  {p(1)=3, D(1)=5, D(k+1)=D(k)+i(k+1), with i(k+1)=4p(k+1)=i(k)+8}
  {Lucas test is done for the first p with mp_jacobi(p^2-4,a)=-1.}
  {This sequence is from Wei Dai / Marcel Martin, see references.}

  b := 3;
  d := 5;
  i := 8;
  repeat
    if mp_jacobi_lm(d,a)=-1 then break;
    if b=bsqr then begin
      {Avoid 'hopeless' loop and test whether a is a perfect square. This }
      {is delayed because we are searching for jacobi(d,a)=-1; if there is}
      {such a d then a is not square, and a square root will be computed  }
      if mp_is_square(a) then goto leave;
    end;
    inc(i,8);
    inc(d,i);
    inc(b,2);
    {assert(d+4=sqr(longint(b)), 'd = b^2 - 4 in mp_is_slpsp_alt');}
    if mp_error<>MP_OKAY then goto leave;
  until b>=bmax;

  {exit if no proper p=b found}
  if b>=bmax then goto leave;

  {p=b>0 should be MUCH less than a but ...}
  mp_set_int(p,b);
  if mp_cmp(a,p)=MP_LT then mp_mod(p,a,p);

  {calculate t,s with a+1 = 2^i*t}
  mp_add_d(a,1,t);
  mp_makeodd(t,t,i);

  {Setup Barrett reduction for a}
  mp_reduce_setup(mu, a);

  {Calculate Lucas sequence}
  s_mp_lucasvmod1(p,a,t,mu,v);

  {t=n-2 for the next compares}
  mp_sub_d(a,2,t);

  if (mp_cmp_d(v,2)=MP_EQ) or (mp_cmp(v,t)=MP_EQ) then begin
    {a is lpsp}
    mp_is_slpsp_alt := true;
    goto leave;
  end;

  for d:=1 to i-1 do begin
    mp_sqr(v,v);      mp_reduce(v,a,mu);
    mp_sub_d(v,2,v);  if v.sign=MP_NEG then mp_add(v,a,v);
    if mp_cmp_d(v,2)=MP_EQ then goto leave;
    if mp_cmp(v,t)=MP_EQ then begin
      {a is lpsp}
      mp_is_slpsp_alt := true;
      goto leave;
    end;
  end;

leave:
  if mp_error<>MP_OKAY then mp_is_slpsp_alt := false;
  mp_clear4(mu,p,t,v);

end;


{---------------------------------------------------------------------------}
function mp_is_slpsp(const a: mp_int): boolean;
  {-Lucas (strong) probable prime test for n using Selfridge A method}
begin
  mp_is_slpsp := s_mp_is_lpsp(a, true);
end;


{---------------------------------------------------------------------------}
function mp_is_lpsp(const a: mp_int): boolean;
  {-Lucas probable prime test for n using Selfridge A method}
begin
  mp_is_lpsp := s_mp_is_lpsp(a, false);
end;


{---------------------------------------------------------------------------}
function s_mp_is_lpsp(const n: mp_int; strong: boolean): boolean;
  {-Lucas (strong) probable prime test for n using Selfridge A method}
var
  q,d,j,s: longint;
  u,v,k,qk,t,mu: mp_int;
const
  dmax = 200000;
begin
  {Ref: - Dana Jacobsen, Pseudoprime Statistics, Tables, and Data, 3 August 2015                  }
  {       http://ntheory.org/pseudoprimes.html and https://github.com/danaj/Math-Prime-Util-GMP   }
  {     - https://en.wikipedia.org/wiki/Lucas_pseudoprime#Implementing_a_Lucas_probable_prime_test}

  {Systematically verified against the tables up to 100000000,}
  {see the test program t_lpsp.pas in the archive supptest.zip}

  {initialize function result with false}
  s_mp_is_lpsp := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_is_lpsp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Handle cases n<=2}
  j := mp_cmp_d(n,2);
  if j=MP_EQ then begin
    s_mp_is_lpsp := true;
    exit;
  end
  else if (j=MP_LT) or mp_iseven(n) then exit;

  {Select D using Selfridge A method}
  d := 5;
  while d < dmax do begin
    j := mp_jacobi_lm(d,n);
    if mp_error<>MP_OKAY then exit;
    if j<1 then break;
    if (d=25) and mp_is_square(n) then exit;
    if d < 0 then d := abs(d) + 2
    else d := -(d+2);
  end;

  if j=0 then begin
    {Jacobi=0 iff gcd(n,d) > 1}
    exit;
  end;
  q := (1 - d) div 4;

  mp_init6(u,v,k,qk,t,mu);
  if mp_error<>MP_OKAY then exit;

  mp_add_d(n,1,k);
  if strong then begin
    {n+1=2^s*k}
    mp_makeodd(k,k,s);
  end;

  j := mp_bitsize(k);
  mp_set1(u);
  mp_set1(v);
  mp_set_int(qk, q);

  {Setup Barrett reduction for n. Note that Barret is actually}
  {SLOWER for small numbers with less than about 100 bits}
  mp_reduce_setup(mu, n);

  while j > 1 do begin
    {U(2k) = U(k) * V(k)}
    mp_mul(u,v,u);
    mp_reduce(u,n,mu);
    {V(2k) = V(k)^2 - 2*Q(k)}
    mp_sqr(v, v);
    mp_mul_2k(qk,1,t);
    mp_sub(v, t, v);
    mp_reduce(v,n,mu);
    {Q(2k) = Q(k)^2}
    mp_sqr(qk, qk);
    mp_reduce(qk,n,mu);
    dec(j);
    if mp_isbit(k,j-1) then begin
      {save D*U(2k) in t}
      mp_mul_int(U, D, t);
      {Note that we have P=1}
      {U2k+1 = (P*U(2k) + V(2k))/2}
      mp_add(u, v, u);
      if mp_isodd(u) then mp_add(u, n, u);
      mp_shr1(u);
      {V(2k+1) = (D*U(2k) + P*V(2k))/2}
      mp_add(v, t, v);
      if mp_isodd(v) then mp_add(v, n, v);
      mp_shr1(v);
      {Q(2k+1)=Q(2k)*Q}
      mp_mul_int(qk, q, qk);
    end;
    mp_reduce(qk,n,mu);
  end;
  mp_reduce(u,n,mu);
  mp_reduce(v,n,mu);

  if mp_error=MP_OKAY then begin
    if strong then begin
      if u.used = 0 then s_mp_is_lpsp := true
      else begin
        while s > 0 do begin
          if v.used=0 then begin
            s_mp_is_lpsp := true;
            break;
          end;
          dec(s);
          if s<>0 then begin
            {V(2k) = V(k)^2 - 2*Q(k)}
            mp_sqr(v, v);
            mp_mul_2k(qk,1,t);
            mp_sub(v, t, v);
            mp_reduce(v,n,mu);
            {Q(2k) = Q(k)^2}
            mp_sqr(qk, qk);
            mp_reduce(qk,n,mu);
          end;
        end;
      end;
    end
    else begin
      {Simple test U(n+1)=0}
      s_mp_is_lpsp := u.used=0;
    end;
  end;
  mp_clear6(u,v,k,qk,t,mu);
end;


{---------------------------------------------------------------------------}
function mp_is_spsp(const n,a: mp_int): boolean;
  {-Strong probable prime test of n to base a > 1 from HAC p. 139 Alg.4.24  }
  { Sets result to false if definitely composite or true if probably prime.}
  { Randomly the chance of error is <= 1/4 and often very much lower.      }
label
  leave;
var
  n1, y, r: mp_int;
  s,j: longint;
  useb: boolean;
const
  SBMin = 5;
begin
  {init default result}
  mp_is_spsp := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_spsp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_cmp_d(n, 1)<>MP_GT then exit;

  {ensure a > 1}
  if mp_cmp_d(a, 1)<>MP_GT then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_is_spsp: a<=1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init3(n1,y,r);
  if mp_error<>MP_OKAY then exit;

  {get n1 = n - 1}
  mp_sub_d(n,1,n1);

  {set 2^s * r = n1}
  mp_makeodd(n1,r,s);

  {compute y = a^r mod n}
  mp_exptmod(a, r, n, y);

  {if y<>1 and y<>n-1 do}
  if (not mp_is1(y)) and mp_is_ne(y, n1) then begin
    j := 1;
    {Use Barret if s is not very small}
    useb := s>=SBMIN;
    if useb then begin
      {Setup Barrett reduction for n}
      mp_reduce_setup(r, n);
    end;
    while (j <= s-1) and mp_is_ne(y, n1) and (mp_error=MP_OKAY) do begin
      if useb then begin
        mp_sqr(y,y);
        mp_reduce(y,n,r);
      end
      else mp_sqrmod(y, n, y);
      {if y=1 then composite}
      if mp_is1(y) then goto leave;
      inc(j);
    end;
    {if y<>n1 then composite}
    if mp_is_ne(y, n1) then goto leave;
  end;

  {probably prime now}
  mp_is_spsp := true;

leave:
  mp_clear3(n1,y,r);
end;


{---------------------------------------------------------------------------}
function mp_is_spsp_d(const n: mp_int; a: mp_digit): boolean;
  {-Strong probable prime test of n to mp_digit base a > 1 from HAC p. 139 Alg.4.24}
var
  t: mp_int;
begin
  mp_is_spsp_d := false;
  if mp_error=MP_OKAY then begin
    mp_init_set(t,a);
    if mp_error=MP_OKAY then begin
      mp_is_spsp_d := mp_is_spsp(n,t);
      mp_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_is_psp(const n,a: mp_int): boolean;
  {-Probable prime test of n to base a > 1, true if a^(n-1) mod n = 1}
var
  n1, y: mp_int;
begin
  {init default result}
  mp_is_psp := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_psp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_cmp_d(n, 1)<>MP_GT then exit;

  {ensure a > 1}
  if mp_cmp_d(a, 1)<>MP_GT then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_is_psp: a <= 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init2(n1,y);
  if mp_error=MP_OKAY then begin
    mp_sub_d(n,1,n1);
    mp_exptmod(a,n1,n,y);
    mp_is_psp := mp_is1(y);
    mp_clear2(n1,y);
  end;
end;


{---------------------------------------------------------------------------}
function mp_is_psp_d(const n: mp_int; a: mp_digit): boolean;
  {-Probable prime test of n to base a > 1, true if a^(n-1) mod n = 1}
var
  t: mp_int;
begin
  mp_is_psp_d := false;
  if mp_error=MP_OKAY then begin
    mp_init_set(t,a);
    if mp_error=MP_OKAY then begin
      mp_is_psp_d := mp_is_psp(n,t);
      mp_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function s_mp_is_psp2(const n: mp_int): boolean;
  {-Test if n>2 is a base-2 (Fermat) probable prime, no init check}
var
  i: longint;
  a,b: mp_int;
  rho: mp_digit;
begin
  s_mp_is_psp2 := false;
  if mp_iseven(n) or (mp_cmp_d(n,2) = MP_LT) then exit;
  mp_init2(a,b);
  s_mp_sub_d(n,1,b);
  mp_montgomery_setup(n, rho);
  mp_montgomery_calcnorm(a, n);
  mp_shl1(a);
  if mp_is_ge(a,n) then mp_sub(a,n,a);
  for i:= mp_bitsize(b)-2 downto 0 do begin
    if MP_Error<>MP_OKAY then break;
    mp_sqr(a,a);
    mp_montgomery_reduce(a, n, rho);
    if mp_isbit(b,i) then begin
      mp_shl1(a);
      if mp_is_ge(a,n) then mp_sub(a,n,a);
    end;
  end;
  mp_montgomery_reduce(a, n, rho);
  s_mp_is_psp2 := mp_is1(a) and (MP_Error=MP_OKAY);
  mp_clear2(a,b);
end;


{---------------------------------------------------------------------------}
function mp_is_pcarmichael_ex(const a: mp_int; numtests: integer; chkprime: boolean): boolean;
  {-Test if a is a Carmichael number with high probability, perform numtests}
  { tests (100 if numtests < 1) and optionally checks if a is prime.}
var
  n: longint;
  k, kmax: integer;
  r,a1,x: mp_int;
const
  cs = 9*25*49*121*169;
  kdef = 100;
begin
  {Default to not Carmichael}
  mp_is_pcarmichael_ex := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_pcarmichael_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Even numbers or numbers < 561 are not Carmichael}
  if mp_iseven(a) or(mp_cmp_d(a,561) < 0) then exit;

  {u >= 561 but fits into longint}
  if mp_is_longint(a,n) then begin
    mp_is_pcarmichael_ex := is_Carmichael32(n);
    exit;
  end;

  {Check square free for 3,5,7,11,13}
  mp_mod_int(a, cs, n);
  n := gcd32(n,cs);
  if (n mod 9 = 0) or (n mod 25 = 0) or (n mod 121 = 0) or (n mod 169 = 0) then exit;

  {Check if a is prime}
  if chkprime then begin
    if mp_is_pprime(a) then exit;
  end;

  mp_init3(r,a1,x);
  if mp_error=MP_OKAY then begin
    {Basic algorithm by Martin Hopf, see}
    {http://math.stackexchange.com/questions/1726016/check-if-a-number-is-carmichael}
    mp_copy(a,a1);
    mp_dec(a1);
    k := 0;
    if numtests<1 then kmax := kdef else kmax := numtests;

    while (k<kmax) and (mp_error=MP_OKAY) do begin
      inc(k);
      repeat
        mp_random(a1,r);
        mp_add_d(r,2,r);
      until mp_is_lt(r,a1);
      mp_exptmod(r,a1,a,x);
      if mp_is1(x) then continue;
      if mp_gcd1(r,a,x) then begin
        {a is composite but not Carmichael}
        break;
      end;
      {a is composite with x a non-trivial divisor of a}
    end;
    mp_is_pcarmichael_ex := (k=kmax) and (mp_error=MP_OKAY);
    mp_clear3(r,a1,x);
  end;
end;


{---------------------------------------------------------------------------}
function mp_is_pcarmichael(const a: mp_int): boolean;
  {-Test if a is is a Carmichael number with high probability, 100 tests, prime check}
begin
  mp_is_pcarmichael := mp_is_pcarmichael_ex(a,0,true);
end;


{---------------------------------------------------------------------------}
function mp_is_square(const a: mp_int): boolean;
  {-Test if a is square}
begin
  mp_is_square := mp_is_square2(a, nil);
end;


{---------------------------------------------------------------------------}
function mp_is_square2(const a: mp_int; psqrt: pmp_int): boolean;
  {-Test if a is square, return sqrt(a) if a is a square and psqrt<>nil}
const {bit i is cleared if i mod n is a square, i=0..m-1}
  ba_37  : array[0..04] of byte = ($64,$e1,$de,$a1,$09);
  ba_41  : array[0..05] of byte = ($c8,$f8,$4a,$7d,$4c,$00);
  ba_43  : array[0..05] of byte = ($ac,$11,$5c,$7c,$a7,$04);
  ba_47  : array[0..05] of byte = ($20,$ac,$d8,$e4,$ca,$7b);
  ba_128 : array[0..15] of byte = ($ec,$fd,$fc,$fd,$ed,$fd,$fd,$fd,$ec,$fd,$fd,$fd,$ed,$fd,$fd,$fd);
const
  DC_MAXLONG = (30+DIGIT_BIT) div DIGIT_BIT;
var
  s,t: mp_int;
  k,r: longint;
  i: integer;
begin
  {Default to Non-square}
  mp_is_square2 := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or ((psqrt<>nil) and mp_not_init(psqrt^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_square2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Negative numbers are not square}
  if a.sign=MP_NEG then exit;

  {0 and 1 are square}
  if (a.used=0) or ((a.used=1) and (a.pdigits^[0]=1)) then begin
    mp_is_square2 := true;
    if psqrt<>nil then mp_copy(a, psqrt^);
    exit;
  end;

  {First check mod 128 (DIGIT_BIT is at least 8)}
  i := a.pdigits^[0] and 127;
  if ba_128[i shr 3] and (1 shl (i and 7)) <> 0 then exit;
  {82.03% rejection rate}

  {Brute force for remaining 32 bit arguments}
  if (a.used<=DC_MAXLONG) and mp_is_longint(a,k) then begin
    if is_square32ex(k,r) then begin
      mp_is_square2 := true;
      if psqrt<>nil then mp_set_int(psqrt^,r);
    end;
    exit;
  end;

  mp_mod_int(a,2029964475,r);   {31*29*25*23*21*17*11}
  if (1 shl (r mod 21)) and $1A7D6C   <> 0 then exit;  {0.6190}
  if (1 shl (r mod 25)) and $D6B5AC   <> 0 then exit;  {0.5600}
  if (1 shl (r mod 31)) and $6DE2B848 <> 0 then exit;  {0.4839}
  if (1 shl (r mod 29)) and $C2EDD0C  <> 0 then exit;  {0.4828}
  if (1 shl (r mod 23)) and $7ACCA0   <> 0 then exit;  {0.4783}
  if (1 shl (r mod 17)) and $5CE8     <> 0 then exit;  {0.4706}
  if (1 shl (r mod 11)) and $5C4      <> 0 then exit;  {0.4545}
  {99.88% cumulative rejection rate}

  mp_mod_int(a,757266679,r); {13*19*37*41*43*47}
  i := r mod 47; if ba_47[i shr 3] and (1 shl (i and 7)) <> 0 then exit; {0.4894}
  i := r mod 43; if ba_43[i shr 3] and (1 shl (i and 7)) <> 0 then exit; {0.4884}
  i := r mod 41; if ba_41[i shr 3] and (1 shl (i and 7)) <> 0 then exit; {0.4878}
  i := r mod 37; if ba_37[i shr 3] and (1 shl (i and 7)) <> 0 then exit; {0.4865}
  if (1 shl (r mod 19)) and $4F50C <> 0 then exit;  {0.4737}
  if (1 shl (r mod 13)) and $9E4   <> 0 then exit;  {0.4615}
  {99.9976% cumulative rejection rate}

  mp_init2(s,t);
  if mp_error=MP_OKAY then begin
    {Final check: test if a-sqrt(a)^2=0}
    s_mp_sqrtrem(a,s,t);
    if mp_is0(t) then begin
      mp_is_square2 := true;
      if psqrt<>nil then mp_exch(s, psqrt^);
    end;
    mp_clear2(s,t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_lucas(k: longint; var lk: mp_int);
  {-Calculate Lucas number lk=luc(k), luc(-k)=(-1)^k*luc(k)}
var
  f1: mp_int;
const
  smluc: array[0..23] of word = (2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,1364,
                                 2207,3571,5778,9349,15127,24476,39603,64079);
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(lk) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lucas');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  { luc(k) = luc(k-1) + luc(k-2) }
  { luc(0) = 2, luc(1)=1         }
  { luc(k) = fib(k) + 2*fib(k-1) }

  if abs(k)<24 then mp_set_w(lk,smluc[abs(k)])
  else begin
    {create temporary for fib(k-1)}
    mp_init(f1);
    if mp_error<>MP_OKAY then exit;
    {get Fibonacci numbers for abs(k)}
    mp_fib2(abs(k),lk,f1);
    mp_mul_2(f1,f1);
    mp_add(lk,f1,lk);
    mp_clear(f1);
  end;
  {if k negative adjust sign of luc(k)}
  if (k<0) and odd(k) then mp_chs(lk,lk);
end;


{---------------------------------------------------------------------------}
procedure mp_lucas2(k: longint; var lk,l1: mp_int);
  {-Calculate two Lucas numbers lk=luc(k), l1=luc(k-1), k>=0}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(lk) or mp_not_init(l1) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lucas2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if k<0 then begin
      {$ifdef MPC_HaltOnArgCheck}
         {$ifdef MPC_UseExceptions}
           raise MPXRange.Create('mp_lucas2: k<0');
         {$else}
           RunError(MP_RTE_RANGE);
         {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
    if @lk=@l1 then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucas2: @lk=@l1');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  {easy outs for k=0,1}
  if k=0 then begin
    mp_set(lk,2);
    mp_set_short(l1, -1);
    exit;
  end
  else if k=1 then begin
    mp_set1(lk);
    mp_set(l1,2);
    exit;
  end;

  {calculate the pair of Lucas numbers from a Fibonacci pair}
  {   L[k] =   F[k] + 2*F[k-1] }
  { L[k-1] = 2*F[k] -   F[k-1] }

  {initialize temporary for F[k-1]}
  mp_init(t); if mp_error<>MP_OKAY then exit;

  {lk = F[k], t = F[k-1]}
  mp_fib2(k,lk,t);

  {L[k-1] = 2*F[k]-F[k-1]}
  mp_mul_2(lk,l1);
  mp_sub(l1,t,l1);

  {L[k] = F[k] + 2*F[k-1]}
  mp_mul_2(t,t);
  mp_add(lk,t,lk);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_lucasv2p(const p,q: mp_int; k: longint; var vk: mp_int; p2: pmp_int; uk: boolean);
  {-Calculate v[k] of Lucas V sequence for p,q; p^2-4q<>0, k>=0;}
  { if p2<>nil, p2^ will be set to v[k+1] or u[k] if uk is true.}
var
  v0,v1,x,t: mp_int;
  i: longint;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) or mp_not_init(q) or mp_not_init(vk) or ((p2<>nil) and mp_not_init(p2^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lucasv2p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if @vk=p2 then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucasv2p: @vk=p2');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  {v[k+2] = p*v[k+1] - q*v[k], v[0]=2, v[1]=p}
  {p^2-4q<>0 is not tested, no proper Lucas sequence!!}

  {easy out for k<=0, k<0 is handled as k=0}
  if k<=0 then begin
    {v[0] = 2}
    mp_set(vk,2);
    if p2<>nil then begin
      if uk then mp_zero(p2^) {u[0] = 0}
      else mp_copy(p,p2^);    {v[1] = p}
    end;
    exit;
  end;

  mp_init4(v0,v1,x,t);
  if mp_error<>0 then exit;

  {v0 = p}
  mp_copy(p, v0);

  {start with x=q}
  mp_copy(q, x);

  {v1 = p*p - 2q}
  mp_sqr(p,v1);
  mp_mul_2(q,t);
  mp_sub(v1,t,v1);

  {loop is never executed for k=1}
  for i:=bitsize32(k)-2 downto 0 do begin
    if mp_error<>MP_OKAY then break;
    if k and (1 shl i) <> 0 then begin
      {v0 = v0*v1 - x*p}
      mp_mul(v0,v1,v0);
      mp_mul(x,p,t);
      mp_sub(v0,t,v0);
      {v1 = v1*v1 - 2x*q}
      mp_sqr(v1,v1);
      mp_mul(x,q,t);
      mp_mul_2(t,t);
      mp_sub(v1,t,v1);
      if i<>0 then begin
        {x = x*x*q}
        mp_sqr(x,x);
        mp_mul(x,q,x);
      end;
    end
    else begin
      {v1 = v0*v1 - x*p}
      mp_mul(v0,v1,v1);
      mp_mul(x,p,t);
      mp_sub(v1,t,v1);
      {v0 = v0*v0 - 2x}
      mp_sqr(v0,v0);
      mp_mul_2(x,t);
      mp_sub(v0,t,v0);
      if i<>0 then begin
        {x = x*x}
        mp_sqr(x,x);
      end;
    end;
  end;

  if mp_error=MP_OKAY then begin
    {store uk[k] or v[k+1] if requested}
    if p2<>nil then begin
      if uk then begin
        if k=1 then begin
          {u[1] = 1}
          mp_set(p2^,1);
        end
        else begin
          {u[k] = (2*v[k+1] - p*v[k]) div (p^2-4q)}
          mp_shl1(v1);
          mp_mul(v0,p,t);
          mp_sub(v1,t,v1);
          mp_sqr(p,x);
          mp_shl(q,2,t);
          mp_sub(x,t,x);
          mp_div(v1,x,p2^);
        end;
      end
      else mp_exch(v1, p2^);
    end;
    {store v[k]}
    mp_copy(v0,vk);
  end;

  mp_clear4(v0,v1,x,t);
end;


{---------------------------------------------------------------------------}
procedure mp_lucasv(const p,q: mp_int; k: longint; var v: mp_int);
  {-Calculate v[k] of Lucas V sequence for p,q, p^2-4q <>0, k>=0}
begin
  mp_lucasv2p(p,q,k,v,nil,false);
end;


{---------------------------------------------------------------------------}
procedure mp_lucasv2(const p,q: mp_int; k: longint; var v,w: mp_int);
  {-Calculate v=v[k],w=v[k+1] of Lucas V sequence for p,q, p^2-4q<>0, k>=0}
begin
  mp_lucasv2p(p,q,k,v,@w,false);
end;


{---------------------------------------------------------------------------}
procedure mp_lucasuv(const p,q: mp_int; k: longint; var uk,vk: mp_int);
  {-Calculate u[k], v[k] of Lucas sequence for p,q; p^2-4q<>0, k>=0}
begin
  mp_lucasv2p(p,q,k,vk,@uk,true);
end;


{---------------------------------------------------------------------------}
procedure mp_lucasu(const p,q: mp_int; k: longint; var u: mp_int);
  {-Calculate u[k] of Lucas sequence for p,q, p^2-4q <>0, k>=0}
var
  v: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  mp_init(v);
  if mp_error=MP_OKAY then begin
    mp_lucasv2p(p,q,k,v,@u,true);
    mp_clear(v);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_lucasvmod(const p,q,n,k: mp_int; var vk: mp_int);
  {-Calculate v[k] mod n of Lucas V sequence for p,q.}
  { Ranges n>1, 0 <= p,q,k < n (checked if MPC_ArgCheck).}
  { Note: p^2-4q<>0 is not tested, no proper Lucas sequence!!}
var
  i: longint;
  v0,v1,x,t,mu: mp_int;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) or mp_not_init(q) or mp_not_init(k) or mp_not_init(n) or mp_not_init(vk) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lucasvmod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if mp_cmp_d(n,2)=MP_LT then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucasvmod: n<2');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    if (p.sign=MP_NEG) or mp_is_ge(p,n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucasvmod: p<0 or p>=n');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    if (q.sign=MP_NEG) or mp_is_ge(q,n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucasvmod: q<0 or q>=n');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    if (k.sign=MP_NEG) or mp_is_ge(k,n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_lucasvmod: k<0 or k>=n');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  {v[k+2] = p*v[k+1] - q*v[k], v[0]=2, v[1]=p}

  {easy out for k<=0, k<0 is handled as k=0}
  if mp_cmp_d(k,1)=MP_LT then begin
    {v[0] = 2}
    mp_set(vk,2);
    exit;
  end;

  mp_init5(v0,v1,x,t,mu);
  if mp_error<>0 then exit;

  {Setup Barrett reduction}
  mp_reduce_setup(mu, n);

  {Paranoia: Initialize with _mod (not really needed if MPC_ArgCheck)}
  {v0 = p}
  mp_mod(p,n,v0);

  {start with x=q}
  mp_mod(q,n,x);

  {v1 = p*p - 2q}
  mp_sqrmod(p,n,v1);
  mp_mul_2(q,t);
  mp_submod(v1,t,n,v1);

  {loop is never executed for k=1}
  for i:=mp_bitsize(k)-2 downto 0 do begin
    if mp_error<>MP_OKAY then break;
    if mp_isbit(k,i) then begin
      {v0 = v0*v1 - x*p}
      mp_mul(v0,v1,v0);     mp_reduce(v0,n,mu);
      mp_mul(x,p,t);        mp_reduce(t,n,mu);
      mp_sub(v0,t,v0);      while (v0.sign=MP_NEG) and (mp_error=MP_OKAY) do mp_add(v0,n,v0);
      {v1 = v1*v1 - 2x*q}
      mp_sqr(v1,v1);        mp_reduce(v1,n,mu);
      mp_mul(x,q,t);        mp_reduce(t,n,mu);
      mp_mul_2(t,t);
      mp_sub(v1,t,v1);      while (v1.sign=MP_NEG) and (mp_error=MP_OKAY) do mp_add(v1,n,v1);
      if i<>0 then begin
        {x = x*x*q}
        mp_sqr(x,x);        mp_reduce(x,n,mu);
        mp_mul(x,q,x);      mp_reduce(x,n,mu);
      end;
    end
    else begin
      {v1 = v0*v1 - x*p}
      mp_mul(v0,v1,v1);     mp_reduce(v1,n,mu);
      mp_mul(x,p,t);        mp_reduce(t,n,mu);
      mp_sub(v1,t,v1);      while (v1.sign=MP_NEG) and (mp_error=MP_OKAY) do mp_add(v1,n,v1);
      {v0 = v0*v0 - 2x}
      mp_sqr(v0,v0);        mp_reduce(v0,n,mu);
      mp_mul_2(x,t);
      mp_sub(v0,t,v0);      while (v0.sign=MP_NEG) and (mp_error=MP_OKAY) do mp_add(v0,n,v0);
      if i<>0 then begin
        {x = x*x}
        mp_sqr(x,x);        mp_reduce(x,n,mu);
      end;
    end;
  end;
  if mp_error=MP_OKAY then mp_exch(v0,vk); {store v[k]}
  mp_clear5(v0,v1,x,t,mu);
end;


{---------------------------------------------------------------------------}
procedure mp_mersenne(n: longint; var mn: mp_int);
  {-Return nth Mersenne number, mn = 2^n-1, MP_RANGE err for n>MaxMersenne}
begin
  mp_2expt(mn, n);
  mp_dec(mn);
end;


{---------------------------------------------------------------------------}
procedure mp_miller_rabin(const n: mp_int; t: word; var prime: boolean);
  {-Miller-Rabin test of n, security parameter t, from HAC p. 139 Alg.4.24}
  { if t=0, calculate t from bit_count with prob of failure < 2^-96}
label
  leave;
var
  n1, y, r: mp_int;
  s,j: longint;
const
  t0: array[0..14] of byte = (50,28,16,10,7,6,5,4,4,3,3,3,3,3,2);
begin
  {default}
  prime := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_miller_rabin');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {easy case n < 2^31 or even}
  if mp_is_longint(n,s) then begin
    if s>1 then prime := IsPrime32(s);
    exit;
  end;
  if mp_iseven(n) or (n.sign=MP_NEG) then exit;

  mp_init3(n1,r,y);
  if mp_error<>MP_OKAY then exit;

  {get n1 = n - 1}
  mp_sub_d(n,1,n1);

  {calculate r,s with n-1=2^s*r}
  mp_makeodd(n1,r,s);

  if t=0 then begin
    {calculate t from bit_count with prob of failure lower than 2^-96}
    j := mp_bitsize(n) div 128;
    if j>14 then t:=14 else t:= word(j);
    t := t0[t];
  end;

  while t>0 do begin
    {generate a in the range 2<=a<n-1, calculate y = a^r mod n}
    repeat
      mp_rand(y, n.used);
      mp_mod(y,n1,y);
      if mp_error<>MP_OKAY then goto leave;
    until mp_cmp_d(y,1)=MP_GT;
    mp_exptmod(y, r, n, y);
    if mp_error<>MP_OKAY then goto leave;

    {if y<>1 and y<>n-1 do}
    if (not mp_is1(y)) and mp_is_ne(y, n1) then begin
      j := 1;
      while (j <= s-1) and mp_is_ne(y, n1) do begin
        mp_sqrmod(y, n, y);
        {if y=1 then composite}
        if mp_is1(y) or (mp_error<>MP_OKAY) then goto leave;
        inc(j);
      end;
      {if y<>n1 then composite}
      if (mp_is_ne(y, n1)) or (mp_error<>MP_OKAY) then goto leave;
    end;
    dec(t);
  end;
  {probably prime now}
  prime := true;

leave:
  mp_clear3(n1,r,y);
end;


{---------------------------------------------------------------------------}
procedure mp_OddProd(a,b: longint; var p: mp_int);
  {-Calculate p=prod(2*i+1),i=a+1...b;  p=1 if a>=b}
const
  MAXL = 26;
var
  st : array[0..MAXL] of mp_int; {local 'stack', avoids init/clear in product}
  lev: integer;                  {stack level}

  procedure OProd(a,b: longint; var p: mp_int);
    {-Calculate f(a+1)*f(a+2)*..f(b-1)*f(b), f(i)=2*i+1}
  const
    dmin={$ifdef BIT16}8{$else}32{$endif}; {min value of d for recursion}
  var
    d: longint;
  begin
    if mp_error<>MP_OKAY then exit;
    d := b-a;
    if d<dmin then begin
      a := 2*b+1;
      mp_set_int(p,a);
      while d>1 do begin
        dec(a,2);
        mp_mul_int(p,a,p);
        dec(d);
      end;
    end
    else begin
      if lev>=MAXL then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXRange.Create('mp_OddProd: lev out of range');
          {$else}
            RunError(MP_RTE_RANGE);
          {$endif}
        {$else}
          set_mp_error(MP_RANGE);
          exit;
        {$endif}
      end;
      inc(lev);
      d := (b+a) shr 1;
      OProd(a,d,p);
      if mp_error=MP_OKAY then begin
        OProd(d,b,st[lev]);
        mp_mul(p,st[lev],p);
      end;
      dec(lev);
    end;
  end;

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_OddProd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if a<b then begin
    lev := -1;
    mp_init_multi(st);
    OProd(a,b,p);
    mp_clear_multi(st);
  end
  else mp_set1(p);
end;


{---------------------------------------------------------------------------}
procedure s_mp_pell4(const d: mp_int; var x,y: mp_int; var r: integer);
  {-Calculate the smallest non-trivial solution of  x^2 - d*y^2 = r}
  { r in [-4,+4,-1,+1]; returns r=0 if d<2 or d is a square. Uses  }
  { continued fraction expansion of sqrt(d) from Forster [4], ch.25}
  { 'Internal' procedure: Assumes d,x,y are initialized and @x<>@y.}
var
  a, w, q, q0, m: mp_int;
  x0, y0, t:  mp_int;
  sign: integer;
  rd: mp_digit;
begin

  r := 0;
  if (mp_error<>MP_OKAY) or (d.used=0) or (d.sign=MP_NEG) then exit;

  {special cases are d=5,12,13. But 13 fits in the algorithm framework}
  {d=5 must be treated separately: solution is x=1, y=1, r=-4}
  if d.used=1 then begin
    rd := d.pdigits^[0];
    if rd=5 then begin
      mp_set1(x);
      mp_set1(y);
      r := -4;
      exit;
    end
    else if rd=12 then begin
      {Forster's original code does not handle the case d=12 separately}
      {and would return (7,2,1), but the (4,1,4) is a smaller solution!}
      mp_set(x,4);
      mp_set1(y);
      r := 4;
      exit;
    end;
  end;

  mp_init8(a, w, q, q0, m, x0, y0, t);

  {calculate w=floor(sqrt(d)), return 0 if d is square}
  s_mp_sqrtrem(d, w, q);
  if mp_is0(q) then begin
    mp_clear8(a, w, q, q0, m, x0, y0, t);
    exit;
  end;

  {Initialize CF algorithm:
   q0 := 1; m := w;
   x0 := 1; x := m;
   y0 := 0; y := 1;}

  mp_set1(q0);  mp_copy(w, m);
  mp_set1(x0);  mp_copy(w, x);
  mp_zero(y0);  mp_set1(y);
  sign := -1;

  {Note 1: the variables with index 1 from [4] are not needed and }
  {swaps are implemented via mp_exch instead of chain assignments.}

  {Note 2: theoretically it is possible to use only one x or y in }
  {the loop instead of both and finally calculate the second from }
  {x^2 - d*y^2 = r. But as long as there is no faster sqrt routine}
  {in this library this alternative is almost always slower.      }

  {loop while q<>1 and q<>4, rd=0 if q=0 or q>MP_DIGIT_MAX, else rd=q}
  if q.used<>1 then rd:=0 else rd := q.pdigits^[0];
  while (mp_error=MP_OKAY) and (rd<>1) and (rd<>4) do begin
    { a := (m + w) div q;}
    mp_add(m,w,a);   mp_div(a,q,a);
    { m1 := a*q -  m;   q1 := q0 + a*(m - m1); }
    { x1 := a*x + x0;   y1 := a*y + y0;        }
    { m  := m1;  q0 := q;    q := q1;          }
    { x0 := x;    x := x1;  y0 := y;  y := y1; }
    mp_mul(a,q,t);   mp_sub(t,m,t);   mp_sub(m,t,m);  mp_exch(m,t);
    mp_mul(a,t,t);   mp_add(t,q0,t);  mp_exch(q,q0);  mp_exch(q,t);
    mp_mul(a,x,t);   mp_add(t,x0,t);  mp_exch(x,x0);  mp_exch(x,t);
    mp_mul(a,y,t);   mp_add(t,y0,t);  mp_exch(y,y0);  mp_exch(y,t);
    {update sign and get next rd}
    sign := -sign;
    if q.used<>1 then rd:=0 else rd := q.pdigits^[0];
  end;

  {if no error calculate final r, avoid D6+ warnings}
  if mp_error=MP_OKAY then r := sign*integer(rd and $7f);
  mp_clear8(a, w, q, q0, m, x0, y0, t);

end;


{---------------------------------------------------------------------------}
procedure mp_pell(const d: mp_int; var x,y: mp_int);
  {-Calculate the smallest non-trivial solution of the Pell equation}
  { x^2 - d*y^2 = 1; error if d<2, if d is a square, or if @x = @y.}
var
  c: integer;
  t,x2: mp_int;
begin
  {Checks in mp_pell4}
  mp_pell4(d,x,y,c);
  if abs(c)=4 then begin
    {calculate ((x+y*sqrt(d))/2)^3}
    mp_init2(t,x2);
    if mp_error=MP_OKAY then begin
      c := c div 4;
      {here c=-1 or 1}
      {x' = x*(x^2-3c)/2}
      {y' = y*(x^2- c)/2}
      mp_sqr(x,x2);
      {t=x2-3*c}
      if c<0 then mp_add_d(x2,3,t) else mp_sub_d(x2,3,t);
      mp_shr(t,1,t);
      mp_mul(x,t,x);
      {t=x2-c}
      if c<0 then mp_add_d(x2,1,t) else mp_sub_d(x2,1,t);
      mp_shr(t,1,t);
      mp_mul(y,t,y);
      mp_clear2(t,x2);
    end;
  end;
  if (c<0) and (mp_error=MP_OKAY) then begin
    {calculate (x+y*sqrt(d))^2, here is x^2 - d*y^2 = -1}
    {x' = 2x^2+1}
    {y' = 2xy   }
    mp_mul(y,x,y);
    mp_shl1(y);
    mp_sqr(x,x);
    mp_shl1(x);
    mp_inc(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_pell4(const d: mp_int; var x,y: mp_int; var r: integer);
  {-Calculate the smallest non-trivial solution of  x^2 - d*y^2 = r,}
  { r in [-4,+4,-1,+1]. Uses continued fraction expansion of sqrt(d)}
  { from Forster [4]. Error if d<2, if d is a square, or if @x = @y.}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(d) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pell4');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @x=@y then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_pell4: @x=@y');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  {call internal procedure, error if r=0}
  s_mp_pell4(d,x,y,r);
  if (mp_error=MP_OKAY) and (r=0) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_pell4: invalid d');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_powerd(const p,q,d: mp_int; n: longint; var x,y: mp_int);
  {-Calculate x + y*sqrt(d) = (p + q*sqrt(d))^n, n >= 0}
var
  r,s,a,b,c: mp_int;
begin
  {See powerd in CALC - a number theory calculator, K.R.Matthews}
  {available from http://www.numbertheory.org/calc/krm_calc.html}

  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) or mp_not_init(p) or mp_not_init(d) or mp_not_init(x) or mp_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_powerd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n<=1 then begin
    if n<1 then begin
      mp_set1(x);
      mp_zero(y);
    end
    else begin
      mp_copy(p,x);
      mp_copy(q,y);
    end;
    exit;
  end;

  mp_init5(r,s,a,b,c);
  if mp_error=MP_OKAY then begin
    {copy const to local before using x,y}
    mp_copy(p,r);
    mp_copy(q,s);
    mp_copy(d,c);
    {right-left powering B = A^n; A = r + s*w;  B = x + y*w,  w = sqrt(d)}
    mp_set1(x);
    mp_zero(y);
    while n>0 do begin
      if odd(n) then begin
        {B := A*B = (r + s*w)*(x + y*w):  x := r*x + c*s*y, y := r*y + s*x}
        mp_mul(r,x,a);
        mp_mul(s,y,b);
        mp_mul(b,c,b);
        mp_add(a,b,b);
        mp_exch(b,x);
        mp_mul(b,s,b);
        mp_mul(y,r,a);
        mp_add(a,b,y);
      end;
      n := n shr 1;
      if n=0 then break;
      {A := A^2 = (r + s*w)^2:  r := r*r + d*s*s,  s := 2rs}
      mp_sqr(r,a);
      mp_sqr(s,b);
      mp_mul(b,c,b);
      mp_add(a,b,b);
      mp_exch(b,r);
      mp_mul(b,s,s);
      mp_shl1(s);
    end;
  end;

  mp_clear5(r,s,a,b,c);
end;


{---------------------------------------------------------------------------}
function s_mp_rqffu(const d: mp_int; var x,y: mp_int): integer;
  {-Return the fundamental unit e = (u + v*sqrt(d))/2 of the real quadratic field}
  { with discriminant d>4, d mod 4 = 0,1. Result is the norm +/- 1 or 0 if error.}
  { Note: This function uses s_mp_pell4 and is much slower than mp_rqffu but the }
  { discriminant may be large. Beware of very long running times for large d!!   }
var
  r: integer;
begin
  s_mp_rqffu := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(d) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_rqffu');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @x=@y then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mp_rqffu: @x=@y');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if (d.used=0) or (d.sign=MP_NEG) or (d.pdigits^[0] and 3 > 1) then begin
    {d <=0 or d mod 4 <>0,1}
    exit;
  end;
  {See Forster[4], Ch. 25, middle of page 245}
  {call internal procedure, error if r=0}
  s_mp_pell4(d,x,y,r);
  if r<>0 then begin
    if abs(r)=1 then begin
      mp_shl1(x);
      mp_shl1(y);
    end
    else begin
      if abs(r)=4 then r := r div 4
      else r := 0;
    end;
  end;
  s_mp_rqffu := r;
end;


{---------------------------------------------------------------------------}
function mp_rqffu(d: longint; var u,v: mp_int): integer;
  {-Return the fundamental unit e = (u + v*sqrt(d))/2 of the real quadratic field}
  { with discriminant d>4, d mod 4 = 0,1. Result is the norm +/- 1 or 0 if error.}
var
  a,p,q,r,t: longint;
  u1,u2,v1,v2,w: mp_int;
begin
  {Ref: Cohen [24], Algorithm 5.7.2 (Fundamental Unit)}
  mp_rqffu := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(u) or mp_not_init(v) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rqffu');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if @u=@v then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_rqffu: @u=@v');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if (d < 1) or (d and 3 > 1) then begin
    {d <= 0 or d mod 4 <>0,1}
    exit;
  end;
  {Step 1: Initialize}
  r := isqrt32(d);
  if d=sqr(r) then exit;
  if (d and 1)=(r and 1) then p := r else p := r-1;
  mp_init5(u1,u2,v1,v2,w);
  mp_set_int(u1,-p);
  mp_set(u2,2);
  mp_set1(v1);
  mp_zero(v2);
  q := 2;
  {Loop exit if error or if period of continued fraction is found (via break)}
  while mp_error=MP_OKAY do begin
    {Step 2: Euclidian step}
    a := (p+r) div q;
    t := p;
    p := a*q-p;
    if (t=p) and (v2.used>0) then begin
      {Step 4: Even period}
      mp_sqr(v2,w);
      mp_mul_int(w,d,w);
      mp_sqr(u2,u);
      mp_add(u,w,u);
      mp_div_int(u,q,@u,t);
      {assert t=0}
      mp_mul(u2,v2,v);
      mp_shl1(v);
      mp_div_int(v,q,@v,t);
      {assert t=0}
      mp_rqffu := 1;
      break;
    end;

    mp_mul_int(u2,a,w);
    mp_exch(u2,u1);
    mp_add(w,u2,u2);
    mp_mul_int(v2,a,w);
    mp_exch(v2,v1);
    mp_add(w,v2,v2);

    {Step 3: Odd period?}
    t := q;
    q := (d - sqr(p)) div q;
    if (t=q) and (v2.used>0) then begin
      mp_mul(v1,v2,w);
      mp_mul_int(w,d,w);
      mp_mul(u1,u2,u);
      mp_add(u,w,u);
      mp_div_int(u,q,@u,t);
      {assert t=0}
      u.sign := MP_ZPOS;
      mp_mul(u1,v2,w);
      mp_mul(u2,v1,v);
      mp_add(v,w,v);
      mp_div_int(v,q,@v,t);
      v.sign := MP_ZPOS;
      {assert t=0}
      mp_rqffu := -1;
      break;
    end;
  end;
  mp_clear5(u1,u2,v1,v2,w);
end;


{---------------------------------------------------------------------------}
procedure mp_poch(n,k: longint; var a: mp_int);
  {-Return the Pochhammer symbol a = (n)_k = n*(n+1)*...(n+k-1), a=1 if k<1}
  { (n)_k is sometimes called "rising factorial" or "Pochhammer function". }
begin
  mp_product(a,n,n+k-1);
end;


{---------------------------------------------------------------------------}
procedure mp_perm(n,r: longint; var a: mp_int);
  {-Compute a=n!/(n-r)!, the number of permutations of n distinct objects}
  { taken r at a time, n,r >= 0.}
var
  k: longint;
  d: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {Init check for a in called routines}
  {$ifdef MPC_ArgCheck}
    if (n<0) or (r<0) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_perm: n<0 or r<0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  {Improved in 1.30.08: handle r=0,1, n-r=0,1. If n-r<=12 and n<=MaxFact}
  {then the quotient of mp_fact and mp_div_d is faster than mp_product. }
  k := n-r;
  if (k<0) or (n<0) then mp_zero(a)
  else if r<=1 then begin
    if r<=0 then mp_set1(a)
    else mp_set_int(a,n);
  end
  else if k<=1 then mp_fact(n,a)
  else if n<=12 then mp_set_int(a,sm_fact[n] div sm_fact[k])
  else begin
    if (k<=12) and (n>5) and (n<=MaxFact) then begin
      k := sm_fact[k];
      if k<=lv_digit_max then begin
        {use only no mp division is needed}
        mp_fact(n,a);
        if k=2 then mp_shr1(a)
        else mp_div_d(a,k,@a,d);
        exit;
      end;
    end;
    mp_product(a,n-r+1,n);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_primorial(n: longint; var a: mp_int);
  {-Primorial of n;  a = n# = product of primes <= n.}
{$ifndef BIT16}
const
  NMAX=$40000;  {max 1MB}
{$else}
const
  NMAX=16000;   {to partial products for n>=176087}
{$endif}
type
  TPTab = array[1..NMAX] of longint;
var
  m,mmax,pps: longint;
  t : mp_int;
  pp: ^TPTab;
  ctx: TPrimeContext;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_primorial');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_set1(a);
  if n<2 then exit;

  {Initialize FindPrime context and allocate memory for array of primes}
  {Prime sieve would be faster but prime generation time with FindNext }
  {is only about 3..4 percent of the total time for n=10^7}
  FindFirstPrime32(2, ctx);

  if n>=mp_primor_cutoff then begin
    mmax := primepi32(n);
    if mmax > NMAX then mmax := NMAX;
    pps := mmax*sizeof(longint);
    pp := IAlloc(pps);
  end
  else begin
    pp   := nil;
    {Keep some compilers quiet}
    pps  := 0;
    mmax := 0;
  end;

  if pp=nil then begin
    {no prime array available, use iterative multiplication}
    while (ctx.prime<=n) and (mp_error=MP_OKAY) do begin
      mp_mul_int(a,ctx.prime,a);
      FindNextPrime32(ctx);
    end;
    exit;
  end;

  {Here we can use the fast mp_prod_int procedure}
  m := 0;
  mp_init(t);
  while (ctx.prime<=n) and (mp_error=MP_OKAY) do begin
    inc(m);
    pp^[m] := ctx.prime;
    if m=mmax then begin
      {array capacity reached, get partial product and accumulate it in a}
      mp_prod_int(t,pp^,m);
      mp_mul(a,t,a);
      {Reset array index}
      m := 0;
    end;
    FindNextPrime32(ctx);
  end;
  {accumulate remaining primes into a}
  if m>0 then begin
    mp_prod_int(t,pp^,m);
    mp_mul(a,t,a);
  end;

  mp_clear(t);
  freemem(pp,pps);
end;


{---------------------------------------------------------------------------}
procedure mp_provable_prime(bits: longint; var p: mp_int);
  {-Generate a random provable prime p with bitsize bits using Maurer's algorithm}
var
  a,q,r: mp_int;

  {-------------------------------------------------------}
  procedure ProvablePrime(k: longint);
    {-Actual recursive core procedure}
  var
    kr: longint;
    f: mp_digit;
    i: integer;
    success: boolean;
  const
    smax = mp_digit(MP_DIGIT_MAX and $3FFF);
    imax = 32;
  begin
    {Generate provable prime with bitsize k using Maurer's algorithm based}
    {on Pocklington's theorem. Simplified version of Alg. 4.62 from HAC[5]}

    {Although ProvablePrime is recursive, the structure of the algorithm}
    {allows the mp_ints to be allocated only once in the outer procedure}

    {Step 1: Special treatment of small bitsizes}
    if k<=31 then begin
      {Step 1.1: Select a random k bit integer}
      mp_rand_bits(p,k);
      {Step 1.2/3: Return next prime}
      {Note that mp_nextprime is exact for bitsizes <= 31}
      mp_nextprime(p);
      exit;
    end;

    {Omit steps 2,3,4: Use fix values B = smax, r = 0.5}
    {Step 5: q = provable_prime(1+floor(r*k)) }
    ProvablePrime(1 + k div 2);

    {Exchange p and q because result of ProvablePrime is returned in p}
    mp_exch(p,q);

    {Get remaining number of bits to generate}
    kr := k - mp_bitsize(q);

    success := false;
    {Step 8}
    while not success do begin
      {Step 8.1: Select next candidate p}
      {try upto imax times to get bitsize(p)=k}
      for i:=1 to imax do begin
        mp_rand_bits(r,kr);
        mp_shl1(r);
        {note that from here r = 2R, R as in HAC, p = 2*R*q + 1}
        mp_mul(r,q,p);
        mp_inc(p);
        if k=mp_bitsize(p) then break;
      end;
      {Step 8.2.a: Use trial division to cast out non-primes}
      mp_small_factor(p, 3, smax, f);
      if f=0 then begin
        {Step 8.2.b: Select random integer >=2}
        mp_rand_bits_ex(a,k-1,false);
        mp_add_d(a,2,a);
        {In the remaining steps 8.2 Pocklington's theorem is used. We have to}
        {compute b1 = a^(n-1) mod p and b2 = a^r mod p. Since n-1=q*r, first}
        {compute b2 = a' = a^r mod p, then b1 = r' = a'^q mod p}
        mp_exptmod(a,r,p,a);
        mp_exptmod(a,q,p,r);
        {check if b1 = r' = a^(n-1) = 1}
        if mp_is1(r) then begin
          {calculate gcd(a^r-1,p) = gcd(a'-1,p)}
          mp_dec(a);
          if mp_gcd1(a,p,r) then begin
            success := true;
          end;
        end;
      end;
    end;
  end;

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_provable_prime');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init3(a,q,r);
  if mp_error=MP_OKAY then begin
    ProvablePrime(bits);
    mp_clear3(a,q,r);
  end;
end;


{---------------------------------------------------------------------------}
function mp_qnr(const n: mp_int): longint;
  {-Return a small quadratic nonresidue for n > 0, -1 if error or no QNR is}
  { found: the smallest prime with kronecker(p|n) = -1 is returned. If n is}
  { a prime, the result is the least positive QNR, for composite n this may}
  { differ from the least QNR: e.g. mp_qnr(15) returns 7 and not 2.}
const
  sq: array[0..9] of integer = (-1,-1,-1, 2,-1, 2,-1, 3, 3,-1);
label
  done;
var
  j,r,m,s: longint;
  i,i0,imax: integer;
  ctx: TPrimeContext;
  big: boolean;
  y: mp_int;

  {-----------------------------------------------------}
  function kronint(p: longint): integer;
    {-Internal Kronecker}
  const
    tab2: array[0..7] of integer = (0,1,0,-1,0,-1,0,1);
  var
    k: integer;
  begin
    if big then begin
      {n > Maxlongint, use n=2^s*y and (p|n) = (p|2)^s * (p|y) }
      {set k = (p|2)^s, skip if s is even}
      if odd(s) then k := tab2[p and 7] else k := 1;
      {(p|2)^s * (p|y)}
      if k=0 then kronint := 0
      else kronint := k*mp_jacobi_lm(p,y);
    end
    else begin
      {n <= MaxLongint, directly use 32-bit function}
      kronint := kronecker32(p,m);
    end;
  end;

begin
  mp_qnr := -1;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_qnr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n.sign=MP_NEG then exit;

  {Last 16-bit prime to check}
  imax := NumPrimes16;

  if mp_is_longint(n,m) then begin
    if m<=9 then begin
      {special cases}
      mp_qnr := sq[m];
      exit;
    end;
    if is_square32(m) then exit;
    if m < $10000 then imax := Primes16Index(m);
    big := false;
  end
  else begin
    mp_init(y);
    if mp_error<>MP_OKAY then exit;
    big := true;
    mp_makeodd(n,y,s);
    if mp_is1(y) then begin
      {n is a power of 2. For odd s 3 is a QNR: (3|2^s) = (3|2)^s = -1;}
      {if s is even n = sqr(2^(s/2)) and there is no QNR for n.}
      if odd(s) then mp_qnr := 3;
      goto done;
    end;
  end;

  {If n is even, start with prime[2] = 3}
  if mp_iseven(n) then i0:=2 else i0:=1;
  for i:=i0 to imax do begin
    r := Primes16[i];
    j := kronint(r);
    if j=-1 then begin
      mp_qnr := r;
      goto done;
    end;
    if (i=25) and big then begin
      {Result > 100, make sure that n is no square}
      if mp_is_square(n) then goto done;
    end;
    if mp_error<>MP_OKAY then goto done;
  end;

  if mp_bitsize(n) > 16 then begin
    FindFirstPrime32($10000,ctx);
    r := ctx.prime;
    while (r>0) and (mp_error=MP_OKAY) do begin
      j := kronint(r);
      if j=-1 then begin
        mp_qnr := r;
        goto done;
      end;
      FindnextPrime32(ctx);
      r := ctx.prime;
    end;
  end;

done:
  if big then mp_clear(y);
end;


{---------------------------------------------------------------------------}
procedure mp_rand_prime(bitsize: longint; pt: TPrimeType; var p: mp_int);
  {-Generate random (probable BPSW) prime of bitsize > 3, pt: prime type of p}
var
  bl: byte;
  done: boolean;
const
  fmaxp = mp_digit(MP_DIGIT_MAX and $3FF);
  Mask7 = mp_digit(MP_DIGIT_MAX and (not 7));
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rand_prime');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if bitsize<4 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rand_prime: bitsize<4');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if mp_error=MP_OKAY then begin
    if (pt=pt_3mod4) or (pt=pt_safe) then bl:=3 else bl:=1;
    repeat
      mp_rand_bits(p,bitsize);
      if pt=pt_1mod8 then p.pdigits^[0] := (p.pdigits^[0] and Mask7) or 1
      else p.pdigits^[0] := p.pdigits^[0] or bl;
      if pt=pt_safe then begin
        {Use mp_safeprime because generating safe primes}
        {with random bits is much less effective.}
        mp_safeprime(p);
        {For small bitsize there are no (or few) safe primes}
        done := (mp_bitsize(p)=bitsize) or (bitsize < 9);
      end
      else done := mp_is_pprime_ex(p, fmaxp);
    until done or (mp_error<>MP_OKAY);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rnr(const a,m: mp_int; var x,y: mp_int);
  {-Rational number reconstruction: for m>0 calculate x,y with a*x=y mod m,}
  { gcd(x,y)=1, 0<=x,|y|<sqrt(m/2), x<>0. x=y=0 if no (unique) solution exists.}
var
  u2,u3,v2,v3,q,sqrtm2: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(m) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rnr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {m must be positive}
  if s_mp_is_le0(m) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rnr: m <= 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init6(u2,u3,v2,v3,q,sqrtm2);
  if mp_error<>MP_OKAY then exit;

  {The following algorithm for rational reconstruction is from Knuth [3],   }
  {Section 4.5.3, Exercise 51. Quote from the answer:                       }

  {'This discussion proves that the problem can be solved efficiently by    }
  {applying Algorithm 4.5.2X with u = m and v = a, but with the following   }
  {replacement for step X2: "If v3 <= sqrt(m/2), the algorithm terminates.  }
  {The pair (x,y) = (|v2|,v3*sign(v2)) is then the unique solution, provided}
  {that gcd(x,y)=1 and x <= sqrt(m/2); otherwise there is no solution."'    }

  {Knuth's Algorithm 4.5.2X is the extended Euclidean algorithm. Here the   }
  {usage of u1,v2 is suppressed and the main loop is simplified via mp_exch }

  {initialize, (v2,v3) = (1,a), force v3 = a mod m into range 0..m-1}
  mp_set(v2,1);
  if (a.sign=MP_NEG) or (mp_cmp_mag(a,m)>=0) then mp_mod(a,m,v3)
  else mp_copy(a,v3);

  {initialize, (u2,u3) = (0,m)}
  mp_set(u2,0);
  mp_copy(m,u3);

  {calculate sqrt(m/2)}
  mp_shr(u3,1,q);
  mp_sqrt(q,sqrtm2);

  {loop invariants  u3>=0, v3>=0}
  while (mp_error=MP_OKAY) and (mp_cmp_mag(v3,sqrtm2)>0) do begin
    mp_divrem(u3, v3, @q, @u3);
    mp_mul(v2, q, q);
    mp_sub(u2, q, u2);
    mp_exch(u2, v2);
    mp_exch(u3, v3);
  end;

  {check if gcd(x,y)=1 and x < sqrt(m/2)}
  if (mp_error=MP_OKAY) and (mp_cmp_mag(v2,sqrtm2)<=0) and mp_gcd1(v2,v3,q) then begin
    mp_abs(v2,x);
    if v2.sign=MP_NEG then mp_chs(v3,y) else mp_copy(v3,y);
  end
  else begin
    mp_zero(x);
    mp_zero(y);
  end;
  mp_clear6(u2,u3,v2,v3,q,sqrtm2);
end;


{---------------------------------------------------------------------------}
procedure mp_rnr2(const a,m,NN,DD: mp_int; var n,d: mp_int);
  {-Rational number reconstruction: for m,NN,DD > 0 calculate co-prime d,n}
  { with a*d=n mod m, |n|<=N, 0<d<=DD, i.e. a=n/d mod m. Return d=n=0 if }
  { no solution exists. The reconstruction is unique if 2*NN*DD < m.}
var
  u2,u3,v2,v3,q: mp_int;
begin
  { This is a generalization of procedure mp_rnr. For a discussion see e.g.}
  { M. Monagan, Maximal quotient rational reconstruction: an almost optimal}
  { algorithm for rational reconstruction. In the Proceedings of ISSAC 2004}
  { pp. 243-249. Available from http://www.cecm.sfu.ca/CAG/papers/MQRR.pdf }

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(m) or mp_not_init(NN) or mp_not_init(DD) or mp_not_init(d) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rnr2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {m must be positive}
  if s_mp_is_le0(m) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rnr2: m <= 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {NN must be positive}
  if s_mp_is_le0(NN) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rnr2: NN <= 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {DD must be positive}
  if s_mp_is_le0(DD) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rnr2: DD <= 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init5(u2,u3,v2,v3,q);
  if mp_error<>MP_OKAY then exit;

  {initialize, (v2,v3) = (1,a), force v3 = a mod m into range 0..m-1}
  mp_set(v2,1);
  if (a.sign=MP_NEG) or (mp_cmp_mag(a,m)>=0) then mp_mod(a,m,v3)
  else mp_copy(a,v3);

  {initialize, (u2,u3) = (0,m)}
  mp_set(u2,0);
  mp_copy(m,u3);

  {loop invariants:  u3>=0, v3>=0}
  while (mp_error=MP_OKAY) and (mp_cmp_mag(v3,NN)>0) do begin
    mp_divrem(u3, v3, @q, @u3);
    mp_mul(v2, q, q);
    mp_sub(u2, q, u2);
    mp_exch(u2, v2);
    mp_exch(u3, v3);
  end;

  {check if gcd(d,n)=1 and |d| <= DD}
  if (mp_error=MP_OKAY) and (mp_cmp_mag(v2,DD)<=0) and mp_gcd1(v2,v3,q) then begin
    mp_abs(v2,d);
    if v2.sign=MP_NEG then mp_chs(v3,n) else mp_copy(v3,n);
  end
  else begin
    mp_zero(d);
    mp_zero(n);
  end;
  mp_clear5(u2,u3,v2,v3,q);
end;


{---------------------------------------------------------------------------}
procedure mp_sigmak(k,n: longint; var a: mp_int);
  {-Compute the sum of the kth powers of the divisors of n; zero if k<0 or n=0.}
  { Special cases k=0: number of divisors of n; k=1: sum of the divisors of n.}
var
  i,j,x: integer;
  v,w: mp_int;
  FR: TFactors32;
begin
  if mp_error<>MP_OKAY then exit;
  {Arg check in zero / set1}

  n := abs(n);
  if (n=0) or (k<0) then begin
    mp_zero(a);
    exit;
  end;

  mp_set1(a);
  if n=1 then exit;

  PrimeFactor32(n, FR);
  with FR do begin
    mp_init2(v,w);
    if mp_error=MP_OKAY then begin
      {sum([primes[i]^(k*(1+pexpo[i])= - 1] / [primes[i]^k - 1], i=1..pcount)}
      {see e.g. NIST[], 27.3.5/6 <http://dlmf.nist.gov/27.3>}
      for i:=1 to pcount do begin
        x := pexpo[i];
        if k=0 then mp_set_int(v,1+x)
        else begin
          {The above formula is not used in order to avoid overflow}
          {of p^[k(x+1)], instead the geometric sum is calculated. }
          if x=1 then begin
            {p^k + 1}
            mp_set_pow(v, primes[i], k);
            mp_inc(v);
          end
          else if x=2 then begin
            {(p^k)^2 + p^k + 1}
            mp_set_pow(w, primes[i], k);
            mp_sqr(w,v);
            mp_add(w,v,v);
            mp_inc(v);
          end
          else begin
            { (((p^k + 1)*p^k + 1)*...)*p^k + 1 with Horner scheme}
            mp_set_pow(w, primes[i], k);
            mp_add_d(w,1,v);
            for j:=1 to x-1 do begin
              mp_mul(v,w,v);
              mp_inc(v);
            end;
          end;
        end;
        mp_mul(a,v,a);
      end;
      mp_clear2(v,w);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_small_factor(const a: mp_int; f0,fmax: mp_digit; var f: mp_digit);
  {-Compute small digit prime factor or 0, f0..fmax, f will be <= min(fmax,$7FFF)}
var
  i,imax,imin: word;
  r: longint;
  q: mp_digit;
const
  s1=3*5*7*11*13*17*19*23;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_small_factor');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {assume no error}
  if fmax<mp_max_small then imax := word(fmax) else imax := word(mp_max_small);

  if f0=0 then imin := 2
  else begin
    if f0>imax then imin:=imax else imin := word(f0);
  end;

  {easy outs}
  f := 0;
  i := a.used;
  if i<1 then exit;
  q := a.pdigits^[0];
  if (i<2) and (q<2) then exit; {a=0 or 1}
  if imin=2 then begin
    if q and 1 = 0 then begin
      f := 2;
      exit;
    end
    else imin := 3;
  end;

  if imin>3 then begin
    {Skip fast test; note that normal (internal) uses have imin 2 or 3}
    i := imin;
    if i and 1 = 0 then inc(i);
  end
  else begin
    mp_mod_int(a,s1,r); {3*5*7*11*13*17*19*23}
    if r mod 3 = 0 then begin
      f := 3;
      exit;
    end;
    if r mod 5 = 0 then begin
      f := 5;
      exit;
    end;
    if r mod 7 = 0 then begin
      f := 7;
      exit;
    end;
    if r mod 11 = 0 then begin
      f := 11;
      exit;
    end;
    if r mod 13 = 0 then begin
      f := 13;
      exit;
    end;
    if r mod 17 = 0 then begin
      f := 17;
      exit;
    end;
    if r mod 19 = 0 then begin
      f := 19;
      exit;
    end;
    if r mod 23 = 0 then begin
      f := 23;
      exit;
    end;
    i := 29;
  end;
  while i<=imax do begin
    if (pbits16[i shr 4] and pmask16[i and $0F] <> 0) then begin
      {Note: i<=imax<=MP_DIGIT_MAX}
      {$ifdef MP_32BIT}
        mp_mod_int(a, i, longint(q));
      {$else}
        mp_mod_d(a, mp_digit(i), q);
      {$endif}
      if mp_error<>MP_OKAY then exit;
      if q=0 then begin
        {return the factor}
        f := mp_digit(i);
        exit;
      end;
    end;
    inc(i,2);
  end;
end;


{---------------------------------------------------------------------------}
{-------------------  Start of mp_sqrtmod routines  ------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure mp_sqrtmod34(const g,p: mp_int; var z: mp_int);
  {-Calculate z with z*z = g mod p, p prime > 0, p mod 4 = 3; *internal*}
var
  t: mp_int;
begin
  mp_init_copy(t,p);
  if mp_error=MP_OKAY then begin
    {Algorithm from P1363/D8 [25], A.2.5, I.}
    mp_inc(t);
    {t = (p+1)/4}
    mp_shr(t,2,t);
    {g^(1/2) = g^((p+1)/4) mod p}
    mp_exptmod(g,t,p,z);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmod58(const g,p: mp_int; var z: mp_int);
  {-Calculate z with z*z = g mod p, p prime = 5 mod 8; *internal*}
var
  r,s,t: mp_int;
begin
  mp_init3(r,s,t);
  if mp_error=MP_OKAY then begin
    {Algorithm from P1363/D8 [25], A.2.5, II.}
    mp_copy(p,s);
    mp_copy(g,r);           {r = g}
    mp_mul_2(g,t);          {t = 2g}
    mp_sub_d(s,5,s);        {s = (p-5)/8}
    mp_shr(s,3,s);
    mp_exptmod(t,s,p,t);    {t = (2g)^(p-5)/8}
    mp_sqrmod(t,p,s);       {s = 2gt^2}
    mp_mul_2(s,s);
    mp_mulmod(s,r,p,s);
    {z = g*t*(s-1) mod p,   [z = g*gamma*(i-1) in P1363]}
    {if s=0 use s-1=p-1}
    if s.used=0 then mp_sub_d(p,1,s) else mp_dec(s);
    mp_mulmod(r,t,p,r);
    mp_mulmod(r,s,p,z);
    mp_clear3(r,s,t);
  end;
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod18_shanks_intern(const a,p: mp_int; var b,q: mp_int; k: longint): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status}
  { *internal* q and k must be setup with p-1 = 2^k*q}
var
  i,j: longint;
  g,r,x: mp_int;
label
  leave;
begin

  {See S.C. Lindhurst [32], Chap.2.2: Shanks's Algorithm for square roots mod p}

  mp_sqrtmod18_shanks_intern := false;
  mp_init3(g,r,x);
  if mp_error=MP_OKAY then begin
    {find quadratic nonresidue g}
    i := 2;
    while mp_jacobi_lm(i,p)<>-1 do begin
      if mp_error<>MP_OKAY then goto leave;
      inc(i);
      if i=100 then begin
        {test p prime and a quadratic residue before going on}
        if (not mp_is_pprime(p)) or (mp_jacobi(a,p)<>1) then goto leave;
      end;
    end;
    mp_set_int(g,i);
    mp_exptmod(g,q,p,g);   {g = g^q}

    mp_div_2(q,q);         {q = (q-1)/2 = q/2, since q is odd}
    mp_exptmod(a,q,p,x);   {x = a^((q-1)/2}
    mp_mulmod(a,x,p,r);    {r = a*x=a^((q+1)/2}
    mp_mulmod(x,r,p,x);    {x = a^q}

    while (not mp_is1(x)) and (mp_error=MP_OKAY) do begin
      {loop invariant: x*r}
      j := 0;
      {find least positive integer with x^(2^j) = 1}
      mp_copy(x,q);
      repeat
        inc(j);
        {j=k should not happen if p is prime, indicate failure}
        if (j=k) or (mp_error<>MP_OKAY) then goto leave;
        mp_sqrmod(q,p,q);
      until mp_is1(q);

      {t := g^(2^(k-j-1)}
      for i:=2 to k-j do mp_sqrmod(g,p,g);

      {update; r=r*g, g=t^2, x=g*x, k=j}
      mp_mulmod(r,g,p,r);
      mp_sqrmod(g,p,g);
      mp_mulmod(x,g,p,x);
      k := j;
    end;
    mp_exch(r,b);
    mp_sqrtmod18_shanks_intern := mp_error=MP_OKAY;
leave:
    mp_clear3(g,r,x);
  end;
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod18_shanks(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status; *internal*}
var
  k: longint;
  q: mp_int;
begin
  mp_sqrtmod18_shanks := false;
  mp_init(q);
  if mp_error=MP_OKAY then begin
    mp_sub_d(p,1,q);
    mp_makeodd(q,q,k);     {p-1 = 2^k*q}
    mp_sqrtmod18_shanks := mp_sqrtmod18_shanks_intern(a,p,b,q,k);
    mp_clear(q);
  end;
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod18_lucas(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 8, return success status}
  { *internal* assumes a,p>1, no arg check}
var
  h,k,vk: mp_int;
  i: longint;
  j: integer;
begin
  {sqrt(a) mod p via Lucas chain,  Crandall/Pomerance: [10], Ex.2.31}
  mp_sqrtmod18_lucas := false;
  mp_init3(h,k,vk);
  if mp_error=MP_OKAY then begin
    {calculate h with jacobi(-4a+h^2,p)=-1}
    mp_mul_d(a,4,k);
    for i:=1 to $3FFF do begin
      mp_set_int(vk,sqr(i));
      mp_sub(vk,k,vk);
      j := mp_jacobi(vk,p);
      if j=-1 then begin
        mp_set_int(h,i);
        break;
      end
      else if i=100 then begin
        {test p prime and a quadratic residue before going on}
        {break loop with j=1 indicating failure}
        if (mp_jacobi(a,p)<>1) or (not mp_is_pprime(p)) then break;
      end;
    end;
    if j=-1 then begin
      {calculate Lucas chain for k=(p+1)/2}
      mp_add_d(p,1,k);
      mp_div_2(k,k);
      mp_lucasvmod(h,a,p,k,vk);
      {b = vk/2 mod p}
      if mp_isodd(vk) then mp_add(vk,p,vk);
      {Lucas failed, may be a is a nonresidue or p is not prime.}
      if not mp_iszero(vk) then begin
        mp_div_2(vk,b);
        mp_sqrtmod18_lucas := mp_error=MP_OKAY;
      end;
    end;
    mp_clear3(h,k,vk);
  end;
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod14_mueller(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 1 mod 4, return success status}
  { *internal* assumes a>1,p>1, no argument checks}
var
  PP,mu,u,v: mp_int;
  j: longint;
  t: word;
begin
  mp_sqrtmod14_mueller := false;

  {S. Mller, On the Computation of Square Roots in Finite Fields,}
  {Designs, Codes and Cryptography, 31, pp. 301-312, 2004         }

  {S. Mller, On the Computation of Square Roots and Cube Roots Modulo p. Slides formerly}
  {available from http://math.uwyo.edu/RMMC/2006/course%20material/RMMC-mueller-2.pdf    }

  (** Reformulated Theorem 3.4, resp. Theorem on p.11 of RMMC-mueller-2.pdf}
   **
   ** Suppose p = 1 mod 4 and (Q|p) = 1
   **  - If (Q-4|p) = -1, let P = Q -2. Then V[(p - 1)/4](P, 1) is a
   **    square root of Q.
   **  - If (Q-4|p) = 1, let P = Q*t^2 - 2, where t is such that
   **    (Q*t^2 - 4 | p) = -1.
   **    Then V[(p - 1)/4](P, 1)/t is a square root of Q.
   **)

  {Variable correspondence to Mller:  a => Q, PP => P, t => t}
  mp_init4(PP,mu,u,v);
  if mp_error<>MP_OKAY then exit;

  {Setup Barrett reduction for p, used in s_mp_lucasvmod1}
  mp_reduce_setup(mu,p);

  {No special treatment of case t=1, fits perfectly in for loop}
  for t:=1 to 255 do begin
    {search u = u(t) = a*t^2 - 4 with (u|p) = -1}
    mp_mul_w(a,sqr(t),u);
    mp_sub_d(u,4,u);
    j := mp_jacobi(u,p);
    if (mp_error<>MP_OKAY) or (j=-1) then break;
  end;
  if (mp_error=MP_OKAY) and (j=-1) then begin
    {PP = a*t^2 - 2 = u(t) + 2, force PP into range [0..p-1]}
    mp_add_d(u,2,PP);
    if PP.sign=MP_NEG then mp_add(PP,p,PP)
    else mp_reduce(PP, p, mu);
    {u = (p-1)/4}
    mp_sub_d(p,1,u);
    mp_shr(u,2,u);
    {v = LucasV[u](PP,1) mod p}
    s_mp_lucasvmod1(PP,p,u,mu,v);
    {b = v/t}
    if t=1 then mp_exch(b,v)
    else begin
      mp_set(u,t);
      mp_invmod(u,p,u);
      mp_mulmod(u,v,p,b);
    end;
    mp_sqrtmod14_mueller := mp_error=MP_OKAY;
  end;
  mp_clear4(PP,mu,u,v)
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod916_kong(const a,p: mp_int; var b: mp_int): boolean;
  {-Calculate b with b*b = a mod p, p prime = 9 mod 16, return success status}
  { *internal* assumes a,p>1, no arg check}

  { This procedure implements Algorithm 2 from: F. Kong et al.,    }
  { Improved generalized Atkin algorithm for computing square roots}
  { in finite fields, Information Processing Letters 98 (2006) 1-5 }
var
  r,s,t,d: mp_int;
  ld: longint;
const
  JACMAX=1000;
label
  leave;
begin
  mp_sqrtmod916_kong := false;
  mp_init4(r,s,t,d);
  if mp_error<>MP_OKAY then exit;

  {r = (2a)^((p-9)/16)};
  mp_addmod(a,a,p,s);
  mp_sub_d(p,9,t);
  mp_shr(t,4,t);
  mp_exptmod(s,t,p,r);
  {t = 2ar^2}
  mp_sqrmod(r,p,t);
  mp_mulmod(s,t,p,t);
  {s = t^2}
  mp_sqrmod(t,p,s);
  {if s=1 then ...}
  if mp_is1(s) then begin
    {select d with (d|p)=-1}
    ld := 2;
    repeat
      if mp_jacobi_lm(ld,p)=-1 then break;
      inc(ld);
      if ld=JACMAX then goto leave;
    until false;
    mp_set_int(d,ld);
    {t = r*d^((p-9)/8}
    mp_sub_d(p,9,t);
    mp_shr(t,3,t);
    mp_exptmod(d,t,p,t);
    mp_mulmod(t,r,p,t);
    {s = 2*(t*d)^2*a}
    mp_mulmod(t,d,p,s);
    mp_sqrmod(s,p,s);
    mp_mulmod(s,a,p,s);
    mp_addmod(s,s,p,s);
    {b = t*d*a*(s-1)}
    if s.used=0 then mp_sub_d(p,1,s) else mp_dec(s);
    mp_mulmod(s,a,p,s);
    mp_mulmod(s,d,p,s);
    mp_mulmod(s,t,p,b);
  end
  else begin
    {b = r*a*(t-1)}
    if t.used=0 then mp_sub_d(p,1,t) else mp_dec(t);
    mp_mulmod(t,a,p,t);
    mp_mulmod(r,t,p,b);
  end;
  mp_sqrtmod916_kong := true;

leave:
  mp_clear4(r,s,t,d);
end;


{---------------------------------------------------------------------------}
function mp_sqrtmod18(const a,p: mp_int; var b: mp_int): boolean;
  {-Internal wrapper function for p=1 mod 8}
var
  k,bs: longint;
  q: mp_int;
begin
  mp_sqrtmod18 := false;
  if mp_sqrtmethod=1 then mp_sqrtmod18 := mp_sqrtmod18_shanks(a,p,b)
  else if mp_sqrtmethod=2 then mp_sqrtmod18 := mp_sqrtmod18_lucas(a,p,b)
  else if mp_sqrtmethod=3 then mp_sqrtmod18 := mp_sqrtmod14_mueller(a,p,b)
  else begin
    {Auto mode}
    {$ifdef MPC_UseKONG}
      if p.pdigits^[0] and 15 = 9 then begin
        mp_sqrtmod18 := mp_sqrtmod916_kong(a,p,b);
        exit;
      end;
    {$endif}
    mp_init(q);
    if mp_error=MP_OKAY then begin
      mp_sub_d(p,1,q);
      mp_makeodd(q,q,k);     {p-1 = 2^k*q}
      bs := mp_bitsize(p);

      {V0.7.08: M.Martin [8] uses Lucas if (bs>128) and (k > bs div 8),}
      {PARI/GP [12] uses Cipolla if k(k-1) > 8*bs + 20, cf. Theorem 3.3}
      {in "Square Roots Modulo p" by Gonzalo Tornar¡a, available from}
      {http://www.cmat.edu.uy/~tornaria/pub/Tornaria-2002.pdf}

      {We use the PARI/Tornar¡a test to select between Shanks and Mueller.}
      {Since Mueller is almost always faster than Lucas/Cipolla this is }
      {quite conservative.}

      if k>2+trunc(sqrt(8.0*bs+20.0)) then mp_sqrtmod18 := mp_sqrtmod14_mueller(a,p,b)
      else mp_sqrtmod18 := mp_sqrtmod18_shanks_intern(a,p,b,q,k);
      mp_clear(q);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmod_ex(const a,p: mp_int; ChkJ: boolean; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p, p prime.}
  { If ChkJ=true then check Jacobi symbol (a|p)=1, Err=-1 if check fails.}
  { Err=1 if failure for p=1 mod 8, Err=2 if p even and not 2.}
begin
  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrtmod_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_is0(p) then begin
    Err:=2;
    exit;
  end;
  if s_mp_mod_is0(a,p) then begin
    mp_zero(b);
    exit;
  end;
  if mp_isodd(p) then begin
    if mp_is1(p) then begin
      Err := 1;
      exit;
    end;
    if mp_is0(a) or mp_is1(a) then begin
      {if a=0 or 1, set b=a}
      mp_copy(a,b);
      exit;
    end
    else begin
      if (not ChkJ) or (mp_jacobi(a,p)=1) then begin
        {selected algorithm depends on p mod 8}
        case p.pdigits^[0] and 7 of
            1: if not mp_sqrtmod18(a,p,b) then Err := 1;
            5: mp_sqrtmod58(a,p,b);
          else mp_sqrtmod34(a,p,b);  {cases 3,7; p=3 mod 4}
        end;
      end
      else Err:=-1;
    end;
  end
  else begin
    {even p, set b=a mod 2 if p=2, else error}
    if mp_cmp_d(p,2)=MP_EQ then mp_mod(a,p,b)
    else Err:=2;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmod(const a,p: mp_int; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p, p prime, with Jacobi check.}
  { Err=-1 if (a|p)<>1, Err=1 if failure for p=1 mod 8, Err=2 if p even and <>2}
begin
  mp_sqrtmod_ex(a,p,true,b,Err);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sqrtmod2k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate unique square root b of an odd integer a with b*b = a mod 2^k and }
  { 0<b<2^k/4 for k>3. Err=1 if a<0; =2 for even a; =3 if a<>1 mod min(2^k,8)}
var
  h,t,w: mp_int;
  j,j1,j2: longint;
  am8: mp_digit;
const
  KThresh = 32;
begin
  Err := 1;
  {Err=1 if a is negative}
  if (mp_error<>MP_OKAY) or (a.sign=MP_NEG) then exit;

  {Err=2 for even a}
  if mp_iseven(a) or (k<=0) then begin
    Err := 2;
    exit;
  end;

  {Err=3 if a<>1 mod min(2^k,8)}
  Err := 3;
  {am8 = a mod 8}
  am8 := a.pdigits^[0] and 7;
  if k<3 then begin
    if ((k=1) and odd(am8)) or ((k=2) and (am8 and 3 = 1)) then begin
      mp_set1(b);
      Err := 0;
    end;
    exit;
  end;
  if am8<>1 then exit;

  {First check some easy cases}
  Err := 0;
  if mp_is1(a) then begin
    mp_set1(b);
    exit;
  end;

  mp_init3(w,h,t);
  s_mp_mod_2k(a,k,b);
  if k >= KThresh then begin
    {Check if a mod 2^k is an ordinary square}
    if mp_is_square2(b,@w) then begin
      mp_exch(w,b);
      mp_clear3(w,h,t);
      exit;
    end;
  end;

  {Algorithm is from P1363/D8 [25], A.2.6 Finding Square Roots Modulo a }
  {power of 2. See also Marcel Martin's ISqrtMod2K in [8] or the section}
  {'Powers of 2' in http://www.johndcook.com/quadratic_congruences.html }

  mp_set1(w);
  mp_set1(h);
  for j:=2 to k-2 do begin
    j1 := j+1;
    if mp_isbit(h,j1) <> mp_isbit(b,j1) then begin
      mp_setbit(w,j);
      {h := (h + w * 2^(j+1)) mod 2^k}
      mp_shl(w,j1,t);
      mp_add(t,h,t);
      s_mp_mod_2k(t,k,h);
      j2 := j+j;
      if j2<k then begin
        {h := h - 2^j2}
        if mp_isbit(h,j2) then begin
          {easy case, just clear the bit}
          mp_clrbit(h,j2);
        end
        else begin
          {subtract 2^j2 and correct if h - 2^j2 < 0}
          mp_2expt(t,j2);
          mp_sub(h,t,h);
          if h.sign=MP_NEG then begin
            mp_setbit(t,k);
            mp_clrbit(t,j2);
            mp_add(h,t,h);
          end;
        end;
      end;
    end;
  end;
  {make sure 0 < w < 2^k/4}
  if mp_isbit(w,k-2) then begin
    mp_2expt(t,k-1);
    mp_sub(t,w,b);
  end
  else mp_exch(w,b);
  mp_clear3(w,h,t);
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmod2k(const a: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate a square root b of an integer a with b*b = a mod 2^k.}
  { Return Err=-1 if there is no solution. For odd a the solutions }
  { are normalized to 0 < b < 2^k/4 for k>3.}
var
  t: mp_int;
  s,h: longint;
  d: mp_digit;
label
  leave;
begin
  if k<=0 then begin
    Err := -1;
    exit;
  end;

  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrtmod2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {easy outs}
  if mp_is0(a) then begin
    mp_zero(b);
    exit;
  end;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  {t := a mod 2^k > 0}
  if a.sign=MP_NEG then mp_mod_2k(a,k,t)
  else mp_copy(a,t);

  {Handle odd a and map error (although all error conditions}
  {should have been handled above)}
  if mp_isodd(t) then begin
    s_mp_sqrtmod2k(t,k,b,Err);
    if Err<>0 then Err := -1;
    goto leave;
  end;

  {If a is even the solutions are from Adler & Coury[26], 4-39. If}
  {solvable, the (one) solution with j=0 is calculated, see below.}

  {Decompose a=2^s*t, a>0, s>0}
  s := mp_cnt_lsb(t);

  if s >= k then begin
    {a = t*2(s-k)*2^k mod 2^k, return b=0}
    mp_zero(b);
    goto leave;
  end;

  Err := -1;
  {if s<k is odd there are no solutions ([26] 4-38)}
  if odd(s) then goto leave;

  mp_shr(t, s, t);
  h := s div 2;
  d := t.pdigits^[0];

  if k=s+1 then begin
    {[26] 4-39,iii: there are 2^h solutions: 1 + j*2^(h+1), j=0..2^h-1.}
    {This is wrong (e.g. for sqrt(16) mod 32), but the actual calculation}
    {in [26] shows that there are 2^h solutions 2^h*(1 + 2*j), j=0..2^h-1.}
    Err := 0;
    mp_2expt(b, h);
  end
  else if k=s+2 then begin
    {[26] 4-39,ii: solvable iff t=1 mod 4; there are 2^(h+1) solutions:}
    {+-2^h + j*2^(h+2), j=0..2^h-1}
    if d and 3 = 1 then begin
      Err := 0;
      mp_2expt(b, h);
    end;
  end
  else if k>s+2 then begin
    {[26] 4-39,i: solvable iff t=1 mod 8; there are 2^(h+2) solutions:}
    {2^h*y + j*2^(k-h), j=0..2^h-1, y one of the 4 solutions of x^2=t mod 2^(k-s)}
    if d and 7 = 1 then begin
      s_mp_sqrtmod2k(t,k-s,b,Err);
      if Err=0 then mp_shl(b,h,b)
      else Err := -1;
    end;
  end;
leave:
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmodp2(const a,p: mp_int; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p^2, p prime.}
  { Alias for mp_sqrtmodpk(,,2,,). Err: error code from mp_sqrtmodpk.}
begin
  mp_sqrtmodpk(a,p,2,b,Err);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sqrtmodpk(const a,p: mp_int; k: longint; red: boolean; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p^k, p odd prime.}
  { Err=-1 if (a|p)<>1, Err=1 if failure for p=1 mod 8, Err=2 if p even.}
  { Err=4 if no inverse mod p^(2^i). If red=true, reduce b mod p^k, otherwise}
  { b may be >= p^k; ex: a=22,p=3,k=3 --> b=34 and 34^2 = 22 mod 27}
var
  p2i,ri,x,z: mp_int;
  i,im: integer;
begin
  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_sqrtmodpk');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_iszero(a) or mp_is1(a) or (k<=0) then begin
    {if a=0 or 1, set b=a}
    if k<=0 then mp_zero(b) else mp_copy(a,b);
    exit;
  end;
  if mp_isodd(p) then begin
    mp_init4(p2i,ri,x,z);
    if mp_error<>MP_OKAY then exit;
    {[10] Alg 2.3.11 (Hensel lifting) with f(x)=a-x^2, f'(x)=-2x}
    {Calls function newr repeatedly until 2^i >= k}
    {r0 = sqrt(a) mod p}
    mp_sqrtmod(a,p,ri,Err);
    if (Err=0) and (k>1) then begin
      if mp_is0(ri) then begin
        {zero solution but a<>0: try to divide out p^2}
        mp_sqr(p,p2i);
        mp_divrem(a,p2i,@x,@z);
        if mp_is0(z) then begin
          s_mp_sqrtmodpk(x,p,k-2,true,ri,Err);
          if Err=0 then mp_mul(ri,p,ri);
        end
        else Err:=4;
      end
      else begin
        {k>1}
        im := bitsize32(k)-1;
        {2^im <= k <= 2^(im+1)}
        if k and pred(k)=0 then begin
          {k a power of two: decrement upper bound}
          dec(im);
          {Now k=2^(im+1) and we are done after im+1 loops}
        end;
        for i:=0 to im do begin
          {calculate p^(2^i): copy from p or square p^(2^(i-1))}
          if i=0 then mp_copy(p,p2i) else mp_sqr(p2i,p2i);
          {z = f'(ri)^-1 mod p^(2^i) = -(2ri)^-1 mod p^(2^i)}
          mp_mul_2(ri,z);
          if not mp_invmodf(z,p2i,z) then begin
            Err := 4;
            break;
          end;
          {x = f(ri)/p^(2^i) = (a-ri^2)/p^(2^i)}
          mp_sqr(ri,x);
          mp_sub(a,x,x);
          mp_div(x,p2i,x);
          {x = -xz mod p^(2^i)) = x/(2ri) mod p^(2^i))}
          mp_mulmod(z,x,p2i,x);
          {r(i+1) = ri + x*p^(2^i)}
          mp_mul(x,p2i,x);
          mp_add(x,ri,ri);
        end;
      end;
    end;
    if Err=0 then begin
      {store last ri as result into b, reduce if requested}
      if red then begin
        mp_expt_int(p,k,x);
        mp_mod(ri,x,b);
      end
      else mp_exch(b,ri);
    end;
    mp_clear4(p2i,ri,x,z);
  end
  else Err:=2;
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtmodpk(const a,p: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate square root b < p^k of a with b*b = a mod p^k, p prime.  }
  { Error codes: if p=2: Err=1 if a<0;  Err=-1 if there is no solution }
  { if p<>2: Err=-1 if (a|p)<>1, Err=1: failure for p=1 mod 8, Err=2 if}
  {          p is even, Err=4: no inverse mod p^(2^i)}
begin
  if mp_cmp_d(p,2)=0 then mp_sqrtmod2k(a,k,b,Err)
  else s_mp_sqrtmodpk(a,p,k,true,b,Err);
end;


(*************************************************************************
The 'quadratic' Hensel lifting used in the standard version of mp_sqrtmodpk
is normally (much) faster than the non-quadratic version listed below, but
it has larger numbers p^(2^i) > p^k if k is not a power of 2.

{---------------------------------------------------------------------------}
procedure mp_sqrtmodpk(const a,p: mp_int; k: longint; var b: mp_int; var Err: integer);
  {-Calculate square root b of a with b*b = a mod p^k, p odd prime.}
  { Err=1 if failure for p=1 mod 8, Err=2 if p even.}
var
  pj,r,x,z: mp_int;
  j: integer;
begin
  err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrtmodpk');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_isodd(p) then begin
    if mp_iszero(a) or mp_is1(a) then begin
      {if a=0 or 1, set b=a}
      mp_copy(a,b);
      exit;
    end
    else begin
      mp_init4(pj,r,x,z);
      if mp_error<>MP_OKAY then exit;
      {Non-quadratic Hensel lifting for f(x)=a-x^2, f'(x)=2x}
      {Example in http://en.wikipedia.org/wiki/Hensel's_lemma}
      mp_sqrtmod(a,p,r,Err);       {r = sqrt(a)}
      if (Err=0) and (k>1) then begin
        mp_mul_2(r,z);
        mp_invmod(z,p,z);          {z = f'(r)^-1 mod p = -(2r)^-1 mod p}
        for j:=1 to k-1 do begin
          if j=1 then mp_copy(p,pj)
          else mp_mul(p,pj,pj);
          mp_sqr(r,x);
          mp_sub(a,x,x);
          mp_div(x,pj,@x,nil);     {x = f(rj)*p^(-j) = (a-rj^2)/p^j)}
          mp_mulmod(z,x,pj,x);     {x = -xz mod p^j) = x/(2r) mod p^j)}
          mp_mul(x,pj,x);
          mp_add(x,r,r);           {r = r + x*p^j}
        end;
      end;
      mp_exch(b,r);
      mp_clear4(pj,r,x,z);
    end;
  end
  else begin
    Err:=2;
  end;
end;
**************************************************************************)


{---------------------------------------------------------------------------}
procedure mp_sqrtmodpq(const a,p,q: mp_int; var x,y: mp_int; var Err: integer);
  {-Calculate square roots +x,-x,+y,-y of a mod (pq); p,q primes}
  { If p=q, x is the root from mp_sqrtmodp2 and y is not computed.}
  { For p<>q: Err=4 if gcd(p,q)<>1, otherwise error code from mp_sqrtmod}
var
  c,d,n,r,s: mp_int;
begin
  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(p) or mp_not_init(q) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrtmodpq');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Handle case p=q: only one solution}
  if mp_is_eq(p,q) then begin
    mp_sqrtmodp2(a,p,x,Err);
    exit;
  end;

  mp_init5(c,d,n,r,s);
  {Ref: Alg. 3.44 from HAC [5], slightly rearranged.}
  {First compute extended gcd, c*p + d*q = gcd(p,q) = n. Error if n<>1}
  {Then compute r = sqrt(a) mod p, s = sqrt(a) mod q}
  {Results x = (rdq+scp) mod pq, y = (rdq-scp) mod pq}
  mp_xgcd(p,q,@c,@d,@n);
  if not mp_is1(n) then Err:=4
  else begin
    {calculate roots mod p and mod q, do nothing if error}
    mp_sqrtmod(a,p,r,Err);
    if (Err=0) and (mp_error=MP_OKAY) then begin
      mp_sqrtmod(a,q,s,Err);
      if (Err=0) and (mp_error=MP_OKAY) then begin
        mp_mul(p,q,n);
        mp_mulmod(d,q,n,d);
        mp_mulmod(c,p,n,c);
        mp_mulmod(r,d,n,r);
        mp_mulmod(s,c,n,s);
        {combine results to get roots x,y. The other roots are -x, -y}
        mp_addmod(r,s,n,x);
        mp_submod(r,s,n,y);
      end;
    end;
  end;
  mp_clear5(c,d,n,r,s);
end;
{---------------------------------------------------------------------------}
{---------------------  End of mp_sqrtmod routines  ------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function mp_squad_mod(const a,b,c,p: mp_int; var x1,x2: mp_int): integer;
  {-Solve a^x + b*x + c = 0 mod p, p odd prime; returns number of roots: 0,1,2.}
  { Primality of p is not checked, if 2a has no inverse, b*x + c = 0 is solved.}
var
  d,s,t: mp_int;
  err: integer;
begin
  mp_squad_mod := 0;
  mp_init3(d,s,t);
  if mp_error<>MP_OKAY then exit;

  mp_shl(a,1,t);
  if mp_invmodf(t,p,t) then begin
    {2a has inverse mod p, compute discriminant d = b^2 - 4ac}
    mp_sqr(b,d);
    mp_mul(a,c,s);
    mp_shl(s,2,s);
    mp_sub(d,s,d);
    mp_mod(d,p,d);
    if mp_is0(d) then begin
      {zero discriminant: x1 = -b * (2a)^(-1) mod p}
      mp_mul(t,b,t);
      s_mp_chs(t);
      mp_mod(t,p,x1);
      mp_squad_mod := 1;
    end
    else begin
      mp_sqrtmod(d,p,s,err);
      if err=0 then begin
        {x1/2 = (-b +/- sqrt(d)) * (2a)^(-1) mod p}
        mp_sub(s,b,x1);
        mp_add(s,b,x2);
        s_mp_chs(x2);
        mp_mulmod(x1,t,p,x1);
        mp_mulmod(x2,t,p,x2);
        mp_squad_mod := 2;
        if mp_is_gt(x1,x2) then mp_exch(x1,x2);
      end;
    end;
  end
  else begin
    {2a has no inverse, solve bx + c = 0}
    if mp_invmodf(b,p,t) then begin
      mp_squad_mod := 1;
      s_mp_chs(t);
      mp_mulmod(t,c,p,x1);
    end
    else begin
      {a=b=0 mod p, no solution if c<>0}
      mp_mod(c,p,t);
      if mp_is0(t) then begin
        mp_squad_mod := 1;
        mp_zero(x1);
      end;
    end;
  end;
  mp_clear3(d,s,t);
end;


{---------------------------------------------------------------------------}
function mp_val(const a: mp_int; r: longint): longint;
  {-Return the valuation of a with respect to r, ie. the largest v with}
  { r^v divides a, v=0 if r=0,1,or -1, v=MaxLongint if a=0}
var
  t: mp_int;
  v: longint;
  i: integer;
begin
  mp_val := 0;
  v := abs(r);
  if (v=0) or (v=1) then exit;
  if mp_is0(a) then begin
    mp_val := MaxLongint;
    exit;
  end;

  if v and (v-1) = 0 then begin
    {r is a power of two}
    i:=0;
    while v and 1 = 0 do begin
      inc(i);
      v := v shr 1;
    end;
    mp_val := mp_cnt_lsb(a) div i;
    exit;
  end;

  mp_init(t);
  if mp_error=MP_OKAY then begin
    mp_valrem(a,r,t,v);
    mp_clear(t);
    mp_val := v;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_valrem_d(const a: mp_int; d: mp_digit; var b: mp_int; var v: longint);
  {-Return valuation v of a with respect to d and remainder b: a = b*d^v}
  { internal iterative version, d no power of two; v=MaxLongint, b=0 if a=0}
var
  g,h: mp_digit;
  j: integer;
  b1: boolean;
  s: mp_int;
begin
  {Internal version used especially for MP_16BIT: Avoid division of}
  {large a by mp_ints with small digit counts > 1, additionally BP7}
  {real mode may produce heap overflow with standard mp_valrem.    }

  mp_copy(a,b);

  if mp_is0(a) then begin
    v := MaxLongint;
    exit;
  end;

  mp_init(s);
  if mp_error<>MP_OKAY then exit;

  {compute g = d^j as the maximum power of d < MP_DIGIT_MAX}
  v := 0;
  h := MP_DIGIT_MAX;
  g := 1;
  j := 0;
  while h>=d do begin
    g := g*d;
    inc(j);
    h := h div d;
  end;

  {Divide out powers of g}
  repeat
    mp_div_d(b,g,@s,h);
    if h=0 then begin
      inc(v,j);
      mp_exch(b,s);
    end;
    b1 := (b.used=1) and (b.pdigits^[0]=1)
  until (h<>0) or b1 or (mp_error<>MP_OKAY);

  {test if more factors of d remain}
  if not b1 and (h mod d = 0) then begin
    repeat
      mp_div_d(b,d,@s,h);
      if h=0 then begin
        inc(v);
        mp_exch(b,s);
      end;
    until (h<>0) or (mp_error<>MP_OKAY);
  end;
  mp_clear(s);
end;


{---------------------------------------------------------------------------}
procedure mp_valrem(const a: mp_int; r: longint; var b: mp_int; var v: longint);
  {-Return valuation v of a with respect to r and remainder b: a = b*r^v}
const
  NRP=22;
var
  rp2: array[0..NRP] of mp_int;
  x,y: mp_int;
  i: integer;
  iter: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_valrem');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  v := abs(r);
  if (r=0) or (v=1) or (a.used=0) then begin
    {special cases r=0,1,or -1: a = a*r^0}
    if a.used=0 then v := MaxLongint
    else v := 0;
    mp_copy(a,b);
    exit;
  end;

  if v and (v-1)=0 then begin
    {r is a power of two}
    i:=0;
    while v and 1 = 0 do begin
      inc(i);
      v := v shr 1;
    end;
    {Here |r|=2^v, get number of least significant bits which are zero}
    {and from this v = the number of blocks with log2(|r|) zero bits. }
    v := mp_cnt_lsb(a) div i;
    {shift away v*log2(|r|) zero bits}
    mp_shr(a, v*i, b);
    {Because we have actually divided by |r|^v, adjust sign if necessary}
    if (r<0) and odd(v) then s_mp_chs(b);
    exit;
  end;

  if v < lv_digit_max then begin
    {$ifdef BIT16}
      iter := (sqr(v) < lv_digit_max) or (a.used > MAXDigits div 2);
    {$else}
      iter := sqr(v) < lv_digit_max;
    {$endif}
    if iter then begin
      s_mp_valrem_d(a, mp_digit(v), b, v);
      if (r<0) and odd(v) then s_mp_chs(b);
      exit;
    end;
  end;

  {Cf. the subquadratic s_mp_radix_astr and GMP [15] function mpz_remove: }
  {Keep on dividing by ascending powers r^(2^i) until a non-zero remainder}
  {is found, then by the stored descending powers giving a zero remainder.}
  mp_init3(x,y,rp2[0]);
  if mp_error<>MP_OKAY then exit;
  mp_copy(a,b);
  mp_set_int(rp2[0],r);
  i := 0;
  while (i<NRP) and (mp_error=MP_OKAY) do begin
    mp_divrem(b,rp2[i],@x,@y);
    if y.used <> 0 then break;
    mp_exch(x,b);
    inc(i);
    mp_init(rp2[i]);
    {Test avoids useless squaring producing zero remainder in next iteration}
    if rp2[i-1].used > 1 + b.used div 2 then break
    else mp_sqr(rp2[i-1], rp2[i]);
  end;

  {v = 2^0 + 2^1 + ... 2^(i-1) = 2^i-1}
  v := pred(longint(1) shl i);
  mp_clear(rp2[i]);

  while (i >= 1) and (mp_error=MP_OKAY) do begin
    dec(i);
    mp_divrem(b,rp2[i],@x,@y);
    if y.used=0 then begin
      inc(v, longint(1) shl i);
      mp_exch(b,x);
    end;
    mp_clear(rp2[i]);
  end;
  mp_clear2(x,y);
end;


{---------------------------------------------------------------------------}
function mp_is_power_modpk(a,n,p: mp_int; k: longint): boolean;
  {-Return true if a is a nth power residue mod p^k;}
  { p prime (not checked); gcd(a,p)=1; n, k > 0.    }
var
  u,v,w: mp_int;
begin
  {Adler/Coury [26], Theorem 6.18 for m=p^k}
  mp_is_power_modpk := false;
  if (mp_error<>0) or (k<1) or (mp_cmp_d(n,1)=MP_LT) then exit;
  mp_init3(u,v,w);
  if mp_error=0 then begin
    if mp_gcd1(a,p,v) then begin
      {v := p^k}
      {u := phi(p^k)=p^k-p^(k-1)}
      mp_expt_int(p,k-1,u);
      mp_mul(u,p,v);
      mp_sub(v,u,u);
      {w := phi(v)/gcd(n, phi(v))}
      mp_gcd(n,u,w);
      mp_div(u,w,w);
      {check a^w mod v = 1}
      mp_exptmod(a,w,v,u);
      mp_is_power_modpk := mp_is1(u);
    end;
    mp_clear3(u,v,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_nroot_modp(const a,n,p: mp_int; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod p; p prime (not checked), gcd(n, p-1)=1 or gcd(n, (p-1)/n)=1}
var
  u,v,m: mp_int;
begin
  Err := 1;
  if mp_error<>0 then exit;
  mp_init3(u,v,m);
  if mp_error=MP_OKAY then begin
    mp_sub_d(p,1,v);
    if mp_invmodf(n,v,u) then begin
      {Euler's theorem a^(p-1)=1 mod p gives x=a^u mod p, u=1/n mod (p-1)}
      mp_exptmod(a,u,p,x);
      Err := 0;
    end
    else begin
      {Reduce n mod phi(n), i.e. assume WLOG 0<n<p-1}
      mp_mod(n,v,m);
      if m.used > 0 then begin
        {Test if p=1 mod n and gcd(n, (p-1)/n)=1}
        mp_divrem(v,m,@v,@u);
        if (mp_error=MP_OKAY) and (u.used=0) then begin
          if mp_invmodf(m,v,u) then begin
            mp_exptmod(a,u,p,v);
            {check if v^n = a}
            mp_exptmod(v,n,p,u);
            mp_submod(u,a,p,u);
            if u.used=0 then begin
              Err := 0;
              mp_copy(v,x);
            end
            else Err := 4; {no n power}
          end
          else Err := 3;   {no inverse of n mod (p-1)/n}
        end
        else Err := 2;     {p <> 1 mod n}
      end;
    end;
    mp_clear3(u,v,m);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_nroot_modpk(const a,n,p: mp_int; k: longint; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod p^k; p prime (not checked), k>0, gcd(a,p)=1 and 1/n mod p^k exists}
var
  u,v: mp_int;
begin
  Err := 1;
  if mp_error<>0 then exit;
  if k<2 then begin
    if k=1 then s_mp_nroot_modp(a,n,p,x,Err)
    else Err := 4;  {k<1}
    exit;
  end;
  mp_init2(u,v);
  if mp_error=0 then begin
    if mp_gcd1(a,p,v) then begin
      {Euler's theorem a^phi(p^k)=1 mod p^k gives }
      {x = a^u mod p^k with u = 1/n mod phi(p^k)) }
      {Compute v = p^k, u = phi(p^k) = p^k-p^(k-1)}
      mp_expt_int(p,k-1,u);
      mp_mul(u,p,v);
      mp_sub(v,u,u);
      if mp_invmodf(n,u,u) then begin
        {here u = 1/n mod phi(p^k)}
        mp_exptmod(a,u,v,x);
        Err := 0;
      end
      else Err := 2;  {gcd(n,phi(p^k))<>1}
    end
    else Err := 3;    {gcd(a,p)<>1}
    mp_clear2(u,v);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_nroot_modpq(const a,n,p,q: mp_int; var x: mp_int; var Err: integer);
  {-Solve x^n = a mod pq; p,q primes (not checked), gcd(n, phi(pq))=1; gcd(a,p)=1 if p=q}
var
  u,v: mp_int;
begin
  Err := 1;
  if mp_error<>0 then exit;
  if mp_is_eq(p,q) then begin
    {Case x^n=a mod p^2}
    s_mp_nroot_modpk(a,n,p,2,x,Err);
    exit;
  end;
  {This is the RSA case, x=a^(1/n) mod phi(pq))}
  mp_init2(u,v);
  if mp_error=0 then begin
    mp_sub_d(p,1,v);
    mp_sub_d(q,1,u);
    mp_mul(v,u,v);
    {here v = phi(pq) = (p-1)*(q-1)}
    if mp_invmodf(n,v,u) then begin
      {pq = (p-1)*(q-1) + p + q - 1}
      mp_add(v,p,v);
      mp_add(v,q,v);
      mp_dec(v);
      {here v=pq; u = 1/n mod phi(pq)}
      mp_exptmod(a,u,v,x);
      Err := 0;
    end
    else Err := 2;  {gcd(n,phi(pq))<>1}
    mp_clear2(u,v);
  end;
end;

{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
procedure s_mp_oddswing(n: longint; var a: mp_int);
  {-Calculate odd part of swing number a = swing(n)}
const
  smallodd: array[0..30] of longint = (
              1,1,1,3,3,15,5,35,35,315,63,693,231,3003,429,6435,6435,109395,
              12155,230945,46189,969969,88179,2028117,676039,16900975,1300075,
              35102025,5014575,145422675,9694845);
{$ifdef BIT16}
const
  FMAXPRIME = $F000 div 4;
{$else}
const
  FMAXPRIME = 105097565;  {=primepi(2^31-1)}
{$endif}
type
  TFTab = array[0..FMAXPrime] of longint;
var
  factab: ^TFTab;
var
  i,p,q,fac,n2,n3,tsize,count,sqrtn,imax: longint;
  useps: boolean;
  sieve: TSieve;
begin

  mp_set1(a);

  if n <= 30 then begin
    if n > 2 then mp_set_int(a, smallodd[n]);
    exit;
  end;

  n2    := n div 2;
  n3    := n div 3;
  count := 0;
  sqrtN := isqrt32(n);

  {BIT16: error if tsize(n) >= 2^16 for n ~ 180000}
  imax := primepi32(n);
  tsize := (imax + 8)*sizeof(longint);
  factab  := mp_getmem(tsize);
  if factab=nil then begin
    set_mp_error(MP_NOTINIT);
    exit;
  end;

  if n <= $FFFF then useps := false
  else begin
    useps := true;
    prime_sieve_init(sieve, 3);
  end;

  for i:=2 to imax do begin
    if useps then p := prime_sieve_next(sieve)
    else p := primes16[i];
    if p <= sqrtN then begin
      fac := 1;
      q := n;
      repeat
        q := q div p;
        if odd(q) then fac := fac*p;
      until q<p;
      if fac>1 then begin
        factab^[count] := fac;
        inc(count);
      end;
    end
    else if p <= n3 then begin
      if odd(n div p) then begin
        factab^[count] := p;
        inc(count);
      end;
    end
    else if p > n2 then begin
      factab^[count] := p;
      inc(count);
    end;
  end;
  if count>0 then mp_prod_int(a,factab^,count);

  if useps then prime_sieve_clear(sieve);
  mp_freemem(pointer(factab),  tsize);
end;


{---------------------------------------------------------------------------}
procedure s_mp_odd_fact(n: longint; var a: mp_int);
  {-Recursive odd part of factorial using oddswing}
var
  s: mp_int;
begin
  if n < 2 then begin
    mp_set1(a);
    exit;
  end;
  s_mp_odd_fact(n div 2, a);
  mp_sqr(a,a);
  if mp_error=MP_OKAY then begin
    mp_init(s);
    if mp_error=MP_OKAY then begin
      s_mp_oddswing(n,s);
      mp_mul(a,s,a);
      mp_clear(s);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_fact_swing(n: longint; var a: mp_int);
  {-Calculate a = factorial(n) using Luschny's prime swing method}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_fact_swing');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n <= 100 then mp_product(a,2,n)
  else begin
    s_mp_odd_fact(n,a);
    mp_shl(a, n - popcount32(n), a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_swing(n: longint; var a: mp_int);
  {-Calculate Luschny swing number a = swing(n) = n!/((n div 2)!)^2}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_swing');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_oddswing(n,a);
  mp_shl(a, popcount32(n), a);
end;


{---------------------------------------------------------------------------}
function mp_digitsum(const a: mp_int; radix: word): longint;
  {-Return the digitsum of a in base radix (2..MAXRadix)}
var
  i,im,bs,rsq,rsum: longint;
  d: mp_digit;
  t: mp_int;
const
  NRP=20; {>log_3(2^31)}
  {$ifdef BIT16}
    BSSMALL = 2000;
  {$else}
    BSSMALL = 500;
  {$endif}
label
  done;
var
  rp2: array[0..NRP] of mp_int;

  {--------------------------------------------------------------------}
  procedure recsum(const x: mp_int; level: integer);
    {-Recursive digitsum}
  var
    q,r: mp_int;
    s: longint;
  begin
    if mp_error<>MP_OKAY then exit;
    if level=0 then begin
      if not mp_is_longint(x,s) or (s>=rsq) then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXRange.Create('mp_digitsum: internal error');
          {$else}
            RunError(MP_RTE_RANGE);
          {$endif}
        {$else}
          set_mp_error(MP_RANGE);
          exit;
        {$endif}
      end;
      while s > 0 do begin
        inc(rsum, s mod radix);
        s := s div radix;
      end;
    end
    else begin
      mp_init2(q,r);
      if mp_error=MP_OKAY then begin
        {Get quotient q and remainder r with x = q*radix^(2*level) + r}
        s_mp_divrem(x,rp2[level],@q,@r);
        dec(level);
        {recsum quotient}
        if q.used>0 then recsum(q,level);
        mp_clear(q);
        {recsum remainder}
        if r.used>0 then recsum(r,level);
        mp_clear(r);
      end;
    end;
  end;

begin
  {Basic algorithm for large numbers is from s_mp_radix_astr}
  {Mathematica:  DigitSum[n_, b_:10] := Total[IntegerDigits[n, b]]}
  mp_digitsum := 0;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_digitsum');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {check range of the radix}
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_digitsum: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {easy out if a=0 or radix=2}
  if a.used=0 then exit;
  if radix=2 then begin
    mp_digitsum := mp_popcount(a);
    exit;
  end;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_abs(a, t);
  bs := mp_bitsize(t);
  rsum := 0;

  if bs <= BSSMALL then begin
    {'small' number, use simple div/mod loop}
    while t.used > 0 do begin
      s_mp_div_d(t,radix,@t, d);
      inc(rsum,d);
    end;
    goto done;
  end;

  mp_init_multi(rp2);
  if mp_error=MP_OKAY then begin
    {Pre-compute power table radix^(2i), also get starting level im,}
    {i.e the smallest value with radix^(2*im) >= sqrt(t)}
    rsq := sqr(radix);
    mp_set(rp2[0], rsq);
    for i:=1 to NRP do begin
      if 2*mp_bitsize(rp2[i-1]) > bs+1 then begin
        im := i-1;
        break;
      end;
      mp_sqr(rp2[i-1], rp2[i]);
      im := i;
    end;
    rsq  := sqr(rsq);
    {do the recursive digitsum}
    recsum(t,im);
    mp_clear_multi(rp2);
  end;

done:
  mp_clear(t);
  mp_digitsum := rsum;
end;


{---------------------------------------------------------------------------}
procedure mp_reverse(const a: mp_int; radix: mp_digit; var b: mp_int);
  {-Calculate b with the reversed digits of a in base radix (2..MAXRadix)}
var
  x: mp_int;
  d,z,q,r,f,block: mp_digit;
  i: integer;
  s: word;
begin
  if mp_error<>MP_OKAY then exit;
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_reverse: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_init(x);
  if mp_error=MP_OKAY then begin
    s := a.sign;
    mp_abs(a,x);
    mp_zero(b);
    block := mp_radpow[radix];
    while (mp_error=MP_OKAY) and (x.used>0) do begin
      {get next block}
      mp_div_d(x,block,@x,d);
      {z will be the block reverse of d}
      z := 0;
      if x.used>0 then begin
        {reverse complete block}
        for i:=1 to mp_radexp[radix] do begin
          q := d div radix;
          r := d - q*radix;
          z := z*radix + r;
          d := q;
        end;
        mp_mul_d(b,block,b);
      end
      else begin
        {last block, reverse until q=0}
        f := 1;
        repeat
          q := d div radix;
          r := d - q*radix;
          f := f*radix;
          z := z*radix + r;
          d := q;
        until q=0;
        mp_mul_d(b,f,b);
      end;
      mp_add_d(b,z,b);
    end;
    if b.used <> 0 then b.sign := s;
    mp_clear(x);
  end;
end;


(*
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
*)

{---------------------------------------------------------------------------}
procedure Calc_MaxFact;
  {-Calculate MaxFact with simple iterations based on Stirling formula}
var
  x,m: double;
  i: integer;
begin
  m := MaxMersenne*ln(2.0);
  x := m/16.0;
  for i:=1 to 8 do x := m/(ln(x)-1.0) - 1.0;
  MaxFact := trunc(x);
end;


begin

  Calc_MaxFact;
  MaxFermat := trunc(ln(MaxMersenne)/ln(2.0));

  mp_sqrtmethod := 0;      {0 = Auto, 1: Shanks, 2: Lucas 3: Mueller}

  {Cutoff point for recursive product in mp_primorial}
  {set rough value. For better results use t_tunepr!}

  {$ifdef BIT16}
    mp_primor_cutoff := 400{1000};
  {$else}
    {$ifdef FPC}
       mp_primor_cutoff := 400{1000};
    {$else}
      {$ifdef MP_32BIT}
        mp_primor_cutoff := 1000{2500};
      {$else}
        mp_primor_cutoff := 6000{16000};
      {$endif}
    {$endif}
  {$endif}

end.

{$ifndef _MP_CONF_INC_}

{$define _MP_CONF_INC_}

(*************************************************************************

 DESCRIPTION : Configuration definitions/options for MPArith library

 REQUIREMENTS: BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     21.03.05  W.Ehrhardt  MPC_ArgCheck,MPC_Assert,MPC_ErrFunc,MPC_HaltOnError,MPC_Diagnostic
 0.11     23.08.05  we          MPC_HaltOnArgCheck
 0.12     27.08.05  we          MPC_UseISAAC, MPC_UseISAAC, MPC_Randomize
 0.13     28.08.05  we          MPC_UseExceptions
 0.14     28.08.05  we          _MP_CONF_INC_, fix FPC RESULT
 0.15     10.01.06  we          $define MPC_UseExceptions if FPC_OBJFPC is defined
                                for FPC VER2_0 (-S2 bug is fixed in V 2.0.2)
 0.16     05.08.06  we          $define MPC_Reduce_2k (use mp_reduce_2k in mp_exptmod)
 0.17     05.08.06  we          $define MPC_Reduce_2k (use mp_reduce_2k in mp_exptmod)
 0.18     15.08.06  we          $define MPC_sqrtmod_fp (enables FP() method in mp_sqrtmod)
 0.19     27.08.06  we          $define MPC_ECM_Primetable (Prime table in mp_ecm_factor)
 0.20     29.12.06  we          $define MPC_USE_Assert (debug mode or compiler supported assert)
 0.21     08.05.07  we          MPC_NOHALT, removed MPC_sqrtmod_fp
 0.22     10.09.07  we          MPC_UseKONG
 0.23     25.01.08  we          MPC_Ln2Ln10Tab
 0.24     21.09.08  we          MPC_UseGCD32
 0.25     08.11.08  we          Do NOT define MPC_ErrFunc as debugging default
 0.26     04.01.09  we          MPC_USE_PRC_30
 0.27     21.06.09  we          MPC_TRACE, MPC_USE_ODS
 0.28     05.07.09  we          Removed MPC_USE_PRC_30
 0.29     01.11.09  we          MPC_Old_EleFun
 0.30     12.02.10  we          MPC_MAXRadix64
 0.31     13.06.10  we          MPC_FPrec30K: MPF_MAX_PREC=30000, default for 16 bit
 0.32     15.01.12  we          MPC_FPrec30K forced for BIT16
 0.33     14.02.12  we          MPC_Ln2Ln10Tab renamed to MPC_E1Ln10Tab
 0.34     25.07.12  we          MPC_FPU_ISQRT, use trunc(sqrt()) for isqrt32
 0.35     06.08.12  we          MPC_SmallSieve, max sieve prime = 1908867043
 0.36     20.09.12  we          removed forced MPC_SmallSieve for BIT16
 0.37     12.12.12  we          MPC_PurePascal (if UINT64 available)
 0.38     23.12.12  we          Removed MPC_Old_EleFun
 0.39     06.05.14  we          MPC_PurePascal $ifdef CPUARM
 0.40     16.07.14  we          MPC_SieveNextPrime: use sieve for mp_next_prime
 0.41     02.03.15  we          FPC3x adjustment
 0.42     03.03.15  we          long-used MPC_ACCEPTBIAS listed here
 0.43     11.08.18  we          MPC_NO_EXECPTIONS
**************************************************************************)

{This file is normally included after STD.INC and 'uses' statements, }
{therefore standard definitions can be used here. BUT NOT in mp_types}
{which has to use sysutils for exception handling}


{.$define MPC_PurePascal}       {Do not use inline asm, assumes HAS_UINT64}

{$ifndef MPC_PurePascal}
  {$ifdef BIT64}
    {$define MPC_PurePascal}
  {$endif}
  {$ifdef CPUARM}
    {$define MPC_PurePascal}
  {$endif}
{$endif}


{$ifdef HAS_INT64}
  {$ifndef MP_16BIT}
    {$define MP_32BIT}          {default 32 bit type mp_digit with int64}
  {$endif}
{$else}
  {$define MP_16BIT}            {default 16 bit type mp_digit without int64}
{$endif}


{$ifdef MP_32BIT}
  {$ifdef BIT16}
    {$undef MP_32BIT}           {NO 32 bit digit types for 16 bit compiler!}
    {$ifndef MP_16BIT}
      {$define MP_16BIT}        {use the 16 bit type mp_digit instead!}
    {$endif}
  {$endif}
{$endif}


{.$define MPC_ArgCheck}         {check arguments: initialized, range etc }
                                {RTE/exception if check failed           }
{$define MPC_Assert}            {use assert for critical parts           }
{.$define MPC_Diagnostic}       {include diagnostic code/functions       }
{$define MPC_ECM_Primetable}    {prime table in mp_ecm_factor (BIT32/64) }
{.$define MPC_ErrFunc}          {mp_error is function not integer var    }
{$define MPC_HaltOnError}       {RTE/exceptions on errors vs. mp_error   }
{$define MPC_HaltOnArgCheck}    {RTE/exceptions on arg check failures    }
{.$define MPC_NOHALT}           {override MPC_HaltOnError/HaltOnArgCheck }
{.$define MPC_NOISAAC}          {override MPC_UseISAAC, no ISAAC         }
{.$define MPC_NO_EXECPTIONS}    {Disable exceptions, use run time errors }
{.$define MPC_Randomize}        {randomize the PRN generator at startup  }
                                {undef if predictable 'random numbers'   }
                                {are needed for debugging or testing     }
{.$define MPC_UseTSC}           {TSC/randomize for PRNG seed, not 64-bit }
{.$define MPC_FPrec30K}         {MPF_MAX_PREC=30000, default for 16 bit  }
{.$define MPC_FPU_ISQRT}        {isqrt32=trunc(sqrt()),default for BIT16 }
{.$define MPC_SmallSieve}       {Max sieve prime = 1908867043            }
{$define MPC_SieveNextPrime}    {use sieve for mp_next_prime             }
{$define MPC_E1Ln10Tab}         {use tables for exp(1) and ln(10)        }
{$define MPC_Reduce_2k}         {use mp_reduce_2k in mp_exptmod          }
{$define MPC_UseKONG}           {use Kong if sqrtmethod=0 and p=9 mod 16 }
{$define MPC_UseGCD32}          {use 32 bit GCD in MP GCD routines       }
{$define MPC_UseISAAC}          {use ISAAC random number generator       }
{.$define MPC_ACCEPTBIAS}       {accept small bias in mp_random_radix    }
{.$define MPC_MAXRadix64}       {use MaxRadix=64 instead of 36           }
{.$define MPC_TRACE}            {emit trace info if mp_verbose > 0       }
{.$define MPC_USE_ODS}          {always use OutputDebugStr for WIN32/64  }
                                {trace even for console mode programs    }

{$ifdef Delphi}
  {$define  MPC_UseExceptions}  {use exception instead of RunError       }
{$endif}

{$ifdef VirtualPascal}
  {$define  MPC_UseExceptions}  {use exception instead of RunError       }
{$endif}

{$ifdef FPC}
  {$ifdef FPC_OBJFPC}
    {$define RESULT}            {fix bug in compiler: -S2 does not work in std.inc}
    {$ifdef FPC2Plus}
      {Real Version 2.0 is buggy, but Version 2.0.2+ works OK}
      {$define MPC_UseExceptions} {use exception instead of RunError}
    {$endif}
  {$endif}
{$endif}

{$ifdef debug}
  {define config options for debug mode}
  {$ifndef MPC_ArgCheck}        {$define MPC_ArgCheck}       {$endif}
  {$ifndef MPC_Assert}          {$define MPC_Assert}         {$endif}
  {$ifndef MPC_HaltOnError}     {$define MPC_HaltOnError}    {$endif}
  {$ifndef MPC_HaltOnArgCheck}  {$define MPC_HaltOnArgCheck} {$endif}
  {$ifndef MPC_Diagnostic}      {$define MPC_Diagnostic}     {$endif}
  {$ifndef MPC_TRACE}           {$define MPC_TRACE}          {$endif}
{$endif}


{$ifdef MPC_PurePascal}
  {$ifdef MPC_UseTSC}           {$undef MPC_UseTSC} {$endif}
{$endif}

{$ifdef MPC_NOISAAC}
  {$ifdef MPC_UseISAAC}         {$undef MPC_UseISAAC} {$endif}
{$endif}

{$ifdef MPC_NOHALT}
  {$ifdef MPC_HaltOnError}      {$undef MPC_HaltOnError}     {$endif}
  {$ifdef MPC_HaltOnArgCheck}   {$undef MPC_HaltOnArgCheck}  {$endif}
{$endif}

{$ifdef MPC_NO_EXECPTIONS}
  {$ifdef MPC_UseExceptions}    {$undef MPC_UseExceptions}   {$endif}
{$endif}

{$ifdef BIT16}
  {$define MPC_FPrec30K}
  {$define MPC_FPU_ISQRT}
  {.$define MPC_SmallSieve}     {Max sieve prime = 1908867043}
  {$ifdef  MPC_ECM_Primetable}  {$undef  MPC_ECM_Primetable} {$endif}
{$endif}

{$ifdef MPC_Assert}
  {$ifdef debug}
    {$define MPC_USE_Assert}
  {$endif}
  {$ifdef HAS_Assert}
    {$define MPC_USE_Assert}
  {$endif}
{$else}
  {$ifdef MPC_USE_Assert} {$undef MPC_USE_Assert} {$endif}
{$endif}

{$endif}

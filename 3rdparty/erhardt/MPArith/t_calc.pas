{Test program for mp_calc unit, (c) W.Ehrhardt 2006-2018}

program t_calc;

{$i std.inc}
{$i mp_conf.inc}

{$x+}  {pchar I/O}
{$i+}  {RTE on I/O error}

{$ifdef BIT16}
{$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}


(*
2007-09-04: line terminator ";", redisplays with "_", script mode
2007-09-09: hex, dec, bin, binomial, var=_, checksum
2007-10-24: measure/display calculation time
2007-10-25: sqrtmod
2008-01-30: and/or/xor
2008-04-15: fix parsing of e.g. 'x mod 10',
            show warning if result has more than $F000 chars
2008-06-15: ispprime
2008-09-25: HexLong, removed mem_util
2008-11-04: Initial Starttimer
2008-12-02: Uses BTypes
2009-01-09: display sign with ; line terminator
2009-02-02: cbrtmod in help display
2010-12-30: !! in help display
2011-21-01: Removed length warning for 32/64 bit
2012-07-15: display MaxBit,MaxFact
2012-08-02: updated help display
2012-09-11: updated help display (add: qnr, prime; remove sqr)
2012-10-05: MEM/PACK (BP7 only)
2017-09-12: Extended help ??
2017-09-19: Radix/oct commands
*)


{
TODO:
 - edit last line

 ? mit val()
}

uses
  BTypes, HRTimer,
  {$ifdef MPC_Diagnostic}
     mp_supp,
  {$endif}
  mp_types, mp_base, mp_numth, mp_calc;


{---------------------------------------------------------------------------}
function HexLong(L: longint): mp_string;
  {-Longint as hex string, LSB first}
var
  i: integer;
  s: string[8];
begin
  s := '';
  for i:=0 to 7 do begin
    s := mp_ucrmap[L and $F] + s;
    L := L shr 4;
  end;
  HexLong := s;
end;


{---------------------------------------------------------------------------}
procedure ShowInfo;
begin
  writeln;
  writeln('Operators:  +  -  *  /  div  mod  ^  !  !!  # (primorial)  % (same as mod)');
  writeln;
  writeln('Functions:  abs(a)           and(a,b)        binomial(a,b)   carmichael(a)');
  writeln('            catalan(a)       cbrtmod(a,b)    digitsum(a,b)   eulerphi(a)');
  writeln('            fermat(a)        fib(a)          gcd(a,b)        invmod(a,b)');
  writeln('            iscarmichael(a)  ispprime(a)     jacobi(a,b)     kronecker(a,b)');
  writeln('            lcm(a,b)         luc(a)          maurer(a)       max(a,b)');
  writeln('            mersenne(a)      min(a,b)        moebius(a)      nextpower(a)');
  writeln('            nextprime(a)     or(a,b)         order(a,b)      perm(a,b)');
  writeln('            poch(a,b)        popcount(a)     prevprime(a)    prime(a)');
  writeln('            primepi(a)       primroot(a)     psp(a,b)        qnr(a)');
  writeln('            random(a)        randprime(a)    root(a,b)       safeprime(a)');
  writeln('            sigma(a,b)       spsp(a,b)       sqrt(a)         sqrtmod(a,b)');
  writeln('            swing(a)         val(a,b)        xor(a,b)');

  writeln;
  writeln('Variables:  x  y  z,  Syntax: Var=expression[;] or Var=_[;]');
  writeln;
  writeln('Other    :  line terminator ";" shows only sign/bitsize/chksum/time of result');
  writeln('            ??[string]<enter>:  extended help about commands string*');
  writeln('            bin, dec, oct, hex: set radix for result display');
  writeln('            radix <n>: set radix for result display');
  writeln('            [r]<string>: enter number in radix r, e.g. [7]42 = 30');
{$ifdef BIT16}
  writeln('            mem/pack : display or pack memory');
{$endif}
  writeln('            "_<enter>" re-displays last result');
{$ifdef CPUARM}
  writeln('            ".<enter>" displays time for last calculation; the ARM version');
  writeln('                       does not use TSC and is not high-precision.');
{$else}
  writeln('            ".<enter>" displays time for last calculation');
{$endif}
  writeln;
end;


const
  MAXH = 47;
const
  hs: array[1..MAXH] of string[72] = (
        'abs(a)           Return the absolute value of a',
        'and(a,b)         Return the bitwise AND of a and b',
        'binomial(a,b)    Calculate the binomial coefficient (a choose b)',
        'carmichael(a)    Return the Carmichael function lambda(a)',
        'catalan(a)       Return the a''th Catalan number',
        'cbrtmod(a,b)     Compute a cube root of a mod b, b prime',
        'digitsum(a,b)    Compute the digit sum of a in base b',
        'eulerphi(a)      Return Euler''s totient function phi(a)',
        'fermat(a)        Return the a''th Fermat number = 2^(2^a)+1',
        'fib(a)           Return the a''th Fibonacci number',
        'gcd(a,b)         Compute the greatest common divisor of a and b',
        'invmod(a,b)      Compute c = a^-1 (mod b)',
        'iscarmichael(a)  Test if a is a Carmichael number',
        'ispprime(a)      Test if a is prime (BPSW probable prime if a>2^32)',
        'jacobi(a,b)      Compute the Jacobi/Legendre symbol (a|b)',
        'kronecker(a,b)   Compute the Kronecker symbol (a|b)',
        'lcm(a,b)         Compute least common multiple of a and b',
        'luc(a)           Calculate the a''th Lucas number',
        'maurer(a)        Generate a random provable prime with a bits',
        'max(a,b)         Return the maximum of a and b',
        'mersenne(a)      Return the a''th Mersenne number = 2^a-1',
        'min(a,b)         Return the minimum of a and b',
        'moebius(a)       Return the Moebius function mu(a)',
        'nextpower(a)     Return the next perfect power >= a',
        'nextprime(a)     Return the next prime >= a',
        'or(a,b)          Return the bitwise OR of a and b',
        'order(a,b)       Return the order of a mod b',
        'perm(a,b)        Compute the number of permutations = a!/(a-b)!',
        'poch(a,b)        Return the Pochhammer symbol (a)_b',
        'popcount(a)      Get population count = number of 1-bits in a',
        'prevprime(a)     Return the previous prime <= a',
        'prime(a)         Return the a''th prime',
        'primepi(a)       Return the number of primes <= a',
        'primroot(a)      Compute the smallest primitive root mod a',
        'psp(a,b)         Probable prime test of a to base b',
        'qnr(a)           Return a small quadratic nonresidue for a',
        'random(a)        Return a pseudo-random with < |a| and >= 0',
        'randprime(a)     Generate random (probable BPSW) prime of bitsize a',
        'root(a,b)        Calculate the b''th root of a = floor(a^(1/b))',
        'safeprime(a)     Compute the next safe prime >= a',
        'sigma(a,b)       Compute the sum of the a''th powers of the divisors of b',
        'spsp(a,b)        Strong probable prime test of a to base b',
        'sqrt(a)          Compute the integer square root floor(sqrt(a)),',
        'sqrtmod(a,b)     Calculate modular square root of a mod b',
        'swing(a)         Calculate Luschny swing number of a',
        'val(a,b)         Return the valuation of a with respect to b',
        'xor(a,b)         Return the bitwise XOR of a and b'
    );

{---------------------------------------------------------------------------}
procedure extended_help(s: bstring);
var
  k: integer;
begin
  for k:=1 to MAXH do begin
    if (s='') or (pos(s,hs[k])=1) then begin
      writeln(hs[k]);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure HelpLine;
begin
  writeln('Type "?<enter>" to get some info about commands, "\q" or "quit" to end.');
  writeln('Type "??[str]<enter>" to get info about commands starting with [str]');
end;


var
  evr: TEval;
  EP,i,n: integer;
  HR: THRTimer;
  ctime: double;
  ac: char8;
  print,script: boolean;
  s: mp_string;
  cmd: string[12];
  ss: string[2];
  radix: word;
  rc: string[2];

const
  qq    : string[2] = '??'; {Keep unicode compilers happy}
  space : string[1] = ' ';

begin

  StartTimer(HR);
  mp_init_eval(evr);

  mp_uppercase := true;
  script := false;
  ctime  := 0.0;
  ac := #0;
  rc := 'D';
  radix := 10;

  for i:=1 to paramcount do begin
    if (paramstr(i)='/s') or (paramstr(i)='/S') then script := true;
  end;

  writeln('T_CALC using MPArith V', MP_VERSION, ' [mp_calc]  (c) W.Ehrhardt 2006-2017');
  writeln('Karatsuba cutoffs:  mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff,
          ',   Toom-3 cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff);
  writeln('Burnikel/Ziegler div cutoff = ',mp_bz_cutoff, ',   MaxBit = ', MP_MAXBIT, ',   MaxFact = ', MaxFact);
{$ifdef debug}
  writeln('MaxPrimorial = ', MaxPrimorial);
{$endif}
  if not script then HelpLine;
  writeln;

  repeat
    if not script then write('[',rc,']:=> ');
    readln(s);
    {s := 'ISCARMICHAEL(561)';}
    {s := 'ISCARMICHAEL(4954039956700380001)';}
    while (s<>'') and (s[1]=' ') do delete(s,1,1);
    while (s<>'') and (s[length(s)]=' ') do delete(s,length(s),1);
    if s='' then continue;

    if s='?' then begin
      ShowInfo;
      continue;
    end;

    if pos(qq,s)=1 then begin
      delete(s,1,2);
      while (s<>'') and (s[1]=' ') do delete(s,1,1);
      i := pos(space,s);
      if i>0 then delete(s,i,length(s));
      extended_help(s);
      continue;
    end;

    if s='.' then begin
      writeln('Time = ', ctime:1:3, ' ms');
      continue;
    end;

    cmd := copy(s,1,12);
    for i:=1 to length(cmd) do cmd[i] := upcase(cmd[i]);

    if (cmd='\Q') or (cmd='Q') or (cmd='QUIT') or (cmd='\@') then break
    {$ifdef VER70}
     else if (cmd='MEM') or (cmd='PACK') then begin
       if cmd='PACK' then begin
         s_mp_shrink(evr.Res);
         s_mp_shrink(evr.X);
         s_mp_shrink(evr.Y);
         s_mp_shrink(evr.Z);
       end;
       writeln('Available memory - total: ', memavail, ',  max: ',maxavail);
       continue;
     end
    {$endif}
    else if cmd='HEX' then begin
      rc := 'H';
      radix := 16;
      if script then writeln(cmd);
      continue;
    end
    else if cmd='DEC' then begin
      rc := 'D';
      radix := 10;
      if script then writeln(cmd);
      continue;
    end
    else if cmd='OCT' then begin
      rc := 'O';
      radix := 8;
      if script then writeln(cmd);
      continue;
    end
    else if cmd='BIN' then begin
      rc := 'B';
      radix := 2;
      if script then writeln(cmd);
      continue;
    end;
    if copy(cmd,1,5)='RADIX' then begin
      if length(cmd)=5 then begin
        writeln('Radix = ', radix);
        continue;
      end;
      delete(cmd,1,5);
      while (cmd<>'') and (cmd[1]=' ') do delete(cmd,1,1);
      {$ifdef D12Plus}
        val(string(cmd),n,i);
      {$else}
        val(cmd,n,i);
      {$endif}
      if (i=0) and (n>1) and (n<=MaxRadix) then begin
        radix := n;
        str(n,rc);
      end
      else writeln('Invalid radix: 2 .. ', Maxradix);
      continue;
    end;

    {Echo input in script mode}
    if script then writeln('[Expr] = ',s);

    {Check for trailing ";", i.e. dont print the result}
    if (s<>'') and (s[length(s)]=';') then begin
      delete(s,length(s),1);
      while (s<>'') and (s[length(s)]=' ') do delete(s,length(s),1);
      print := false;
    end
    else print := true;
    if s='' then continue;

    if s<>'_' then begin
      s := s + #0;
      ac := #0;
      if (length(s)>1) and (upcase(s[1]) in ['X','Y','Z']) then begin
        {Check if we have an assignment to x,y,z or just an expression}
        {starting with x,y,z like x mod 10. Analyse first non-blank char}
        {next while/if are OK because s has trailing #0}
        while (s[2]=' ') and (s[3]=' ') do delete(s,2,1);
        if (s[2]=' ') and (s[3]='=') then delete(s,2,1);
        if s[2]='=' then begin
          ac := upcase(s[1]);
          delete(s,1,2);
          while (s<>'') and (s[1]=' ') do delete(s,1,1);
        end;
      end;
      if s<>'_'#0 then begin
        RestartTimer(HR);
        mp_calculate(@s[1],evr,EP);
        ctime := 1000.0*ReadSeconds(HR);
        if evr.Err=Err_MPERR_Eval then set_mp_error(MP_OKAY);
      end;
    end
    else if evr.Err>0 then continue;

    if evr.Err>0 then begin
      writeln('Error ', evr.Err,', ',mp_calc_errorstr(evr.Err),':  <',copy(s,EP+1,length(s)-EP-1), '>');
      if evr.Err=Err_Unknown_Function then HelpLine;
    end
    else if evr.Err<0 then writeln(#7'Error: ', evr.Err, ' [',mp_calc_errorstr(evr.Err),']')
    else begin
      if ac in ['X', 'Y', 'Z'] then begin
        case ac of
          'X':  mp_copy(evr.Res,evr.X);
          'Y':  mp_copy(evr.Res,evr.Y);
          'Z':  mp_copy(evr.Res,evr.Z);
        end;
        write(ac,' = ');
      end
      else write('Result = ');
      case mp_sign(evr.Res) of
        +1 : ss := '>0';
         0 : ss := '=0';
        else ss := '<0'
      end;
      if print then begin
        {$ifdef BIT16}
          if mp_radix_size(evr.Res, radix) > $F000 then begin
            {writeln(mp_radix_astr(evr.Res, radix));}
            writeln('** to many chars in result **');
            write(' [',ss, ', ',mp_bitsize(evr.Res):6, ' bits,  chksum=$',HexLong(mp_checksum(evr.Res)),
                  ',  time=', ctime:1:3, ' ms]');
          end
          else mp_output_radix(evr.Res, radix)
        {$else}
          mp_output_radix(evr.Res, radix)
        {$endif}
      end
      else begin
        write(' [',ss, ', ',mp_bitsize(evr.Res):6, ' bits,  chksum=$',HexLong(mp_checksum(evr.Res)),
              ',  time=', ctime:1:3, ' ms]');
      end;
    end;
    writeln;
  until false;
  mp_clear_eval(evr);
  {$ifdef MPC_Diagnostic}
     mp_dump_meminfo;
     mp_dump_diagctr;
  {$endif}
end.


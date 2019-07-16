{Simple test for mt19937 unit, we Nov.2008}

program t_rnd_44;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$else}
    crt,
  {$endif}
  mt19937,ministat;

const
  NMAX = MaxLongint;
label
  _break;
var
  mx,sx,x,xmin,xmax: double;
  n: longint;
  stat: TStatX;
  rng: mt19937_ctx;
begin
  writeln('Test program for mt19937 unit, mean/sdev of prng_dword    (c) 2008 W.Ehrhardt');
  writeln('MT19937 selftest: ', mt19937_selftest);
  writeln('Count':12, 'Mean':20, 'Min':20, 'Max':20);
  stat1_init(stat);
  mt19937_init(rng,1);
  xmin := 1E50;
  xmax :=-1E50;
  for n:=1 to NMAX do begin
    x := mt19937_dword(rng);
    stat1_add(stat,x);
    if x<xmin then xmin := x;
    if x>xmax then xmax := x;
    if n and $3FFFF = 0 then begin
      stat1_result(stat,mx,sx);
      if stat.Error<>0 then begin
        writeln('Ministat error: ',stat.Error);
        halt;
      end;
      writeln(n:12, mx:20:2, xmin:20:2, xmax:20:2);
      if keypressed and (readkey=#27) then goto _break;
    end;
  end;
_break:
  stat1_result(stat,mx,sx);
  writeln(stat.Nn:12, mx:20:2, xmin:20:2, xmax:20:2);
end.

{Simple test for ISAAC unit, we Nov.2008}

program t_rnd_34;

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
  isaac,ministat;

const
  NMAX = MaxLongint;
label
  _break;
var
  mx,sx,x,xmin,xmax: double;
  n: longint;
  stat: TStatX;
  rng: isaac_ctx;
begin
  writeln('Test program for ISAAC unit, mean/sdev of prng_dword    (c) 2008 W.Ehrhardt');
  writeln('ISAAC selftest: ', isaac_selftest);
  writeln('Count':12, 'Mean':20, 'Min':20, 'Max':20);
  stat1_init(stat);
  isaac_init(rng,1);
  xmin := 1E50;
  xmax :=-1E50;
  for n:=1 to NMAX do begin
    x := isaac_dword(rng);
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

{Simple test for taus88 unit, we Nov.2008}

program t_rnd_14;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$else}
    crt,
  {$endif}
  taus88,ministat;

const
  NMAX = MaxLongint;
label
  _break;
var
  mx,sx,x,xmin,xmax: double;
  n: longint;
  stat: TStatX;
  rng: taus88_ctx;
begin
  writeln('Test program for taus88 unit, mean/sdev of taus88_dword    (c) 2008 W.Ehrhardt');
  writeln('Taus88 selftest: ', taus88_selftest);
  writeln('Count':12, 'Mean':20, 'Min':20, 'Max':20);
  xmin := 1E50;
  xmax :=-1E50;
  stat1_init(stat);
  taus88_init(rng,1);
  for n:=1 to NMAX do begin
    x := taus88_dword(rng);
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

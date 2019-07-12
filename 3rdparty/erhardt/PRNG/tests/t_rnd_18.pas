{Simple test for taus88 unit, we Dec.2017}

program t_rnd_18;

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
  {$endif}
  taus88;
var
  rng: taus88_ctx;
  i: longint;
  w: word;
  buf: array[0..100] of longint;
const
  NTEST = 100*100000;
begin
  writeln('Test for taus88 unit: function rangew(100)    (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', taus88_selftest);
  taus88_init(rng,1);
  fillchar(buf, sizeof(buf), 0);
  for i:=1 to NTEST do begin
    w := taus88_rangew(rng,100);
    if w>99 then begin
      writeln('Error w=', w, '  for i=',i);
    end
    else inc(buf[w]);
  end;
  for i:=0 to 99 do begin
    write(buf[i]:10);
    if (i mod 7 = 6) or (i=99) then writeln;
  end;
end.

{Simple test for xor4096 unit, we Dec.2017}

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
  xor4096;
var
  rng: xor4096_ctx;
  i: longint;
  w: word;
  buf: array[0..100] of longint;
const
  NTEST = 100*100000;
begin
  writeln('Test for xor4096 unit: function rangew(100)    (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', xor4096_selftest);
  xor4096_init(rng,1);
  fillchar(buf, sizeof(buf), 0);
  for i:=1 to NTEST do begin
    w := xor4096_rangew(rng,100);
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

{Part 2 of regression test for DAMath unit,   (c) 2018  W.Ehrhardt}
{Tests for double-double functions}

unit t_damat2;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


interface

procedure test_ddivrem;
procedure test_fma_d;
procedure test_misc2d;


implementation

uses
  damath, t_damat0;


{---------------------------------------------------------------------------}
procedure test_ddivrem;
var
  i,k: integer;
  failed: longint;
  a,b,x,y,f: double;
begin
  writeln('Function ddivrem');
  failed := 0;

  for i:=-50 to 50 do begin
    if i<>0 then begin
      for k:=1 to 1000 do begin
        a := i;
        b := k;
        ddivrem(a,b,x,y);
        f := b*x + y;
        if abs(a-f) > 0.5*eps_d*abs(a) then begin
          inc(failed);
          writeln(i:4, k:4, a-f:30);
        end;
      end;
    end;
  end;

  if failed>0 then begin
    writeln('*** failed');
    inc(total_failed);
  end
  else writeln(' All tests OK');
  inc(total_cnt);
end;


{---------------------------------------------------------------------------}
procedure test_fma_d;
var
  x,y,f: double;
begin
  write('Test fma_d ... ');
  x := ldexpd(1,-30);
  y := fma_d(1+x,1-x,-1);
  f := -x*x;
  if reldev(y,f) > EPS then begin
    writeln('*** failed: ', reldev(y,f));
    inc(total_failed);
  end
  else  writeln('OK.');
  inc(total_cnt);
end;



{---------------------------------------------------------------------------}
procedure test_misc2d;
var
  a,b,c,d: dbl2;
  fail, cnt, k, n: integer;
const
  EPS = 3e-29;
const
  lnpih:  THexDblW = ($A1BD,$48E7,$50D0,$3FF2);  { 1.14472988584940E+0000}
  lnpil:  THexDblW = ($5088,$AD8D,$ABF2,$3C67);  { 1.02659511627078E-0017}
  exppih: THexDblW = ($933A,$6EB0,$2404,$4037);  { 2.31406926327793E+0001}
  exppil: THexDblW = ($1952,$2DD8,$4C96,$BCD8);  {-1.34887470919957E-0015}

  procedure test1(num: integer; x: dbl2; i: longint);
  begin
    inc(cnt);
    if (x.h<>i) or (x.l<>0.0) then begin
      writeln('Test ', num, ' failed: ',i,' <> ',x.h:26,x.l:20);
    end;
  end;

  procedure test2(num: integer; x: dbl2);
  begin
    inc(cnt);
    if abs(x.h+x.l) > EPS then begin
      inc(fail);
      writeln('Test ', num, ' failed: ',x.h:26,x.l:26);
    end;
  end;

begin
  writeln('Test other dbl2 functions ... ');
  cnt := 0;
  fail := 0;

  ddto2d(1, 1e-70,a);
  ceil2d(a,b);
  test1(1, b, 2);
  floor2d(a,b);
  test1(2, b, 1);
  trunc2d(a,b);
  test1(3, b, 1);

  ddto2d(1, -1e-70,a);
  ceil2d(a,b);
  test1(4, b, 1);
  floor2d(a,b);
  test1(5, b, 0);
  trunc2d(a,b);
  test1(6, b, 0);

  ddto2d(-1, 1e-70,a);
  ceil2d(a,b);
  test1(7, b, 0);
  floor2d(a,b);
  test1(8, b, -1);
  trunc2d(a,b);
  test1(9, b, 0);

  ddto2d(-1, -1e-70,a);
  ceil2d(a,b);
  test1(10, b, -1);
  floor2d(a,b);
  test1(11, b, -2);
  trunc2d(a,b);
  test1(12, b, -1);

  dto2d(1,a);
  exp2d(a,b);
  exp1_2d(c);
  sub2d(c,b,c);
  test2(13,c);

  exp1_2d(a);
  ln2d(a,b);
  add21d(b,-1,c);
  test2(14, c);

  pi2d(a);
  ln2d(a,b);
  hhto2d(lnpih, lnpil, d);
  sub2d(b,d,c);
  test2(15,c);
  exp2d(d,b);
  sub2d(b,a,c);
  test2(16,c);

  pi2d(a);
  exp2d(a,b);
  hhto2d(exppih, exppil, d);
  sub2d(b,d,c);
  test2(17,c);
  ln2d(d,b);
  sub2d(b,a,c);
  test2(18,c);

  dto2d(2,a);
  ln2d(a,b);
  ln2_2d(c);
  sub2d(b,c,c);
  test2(19,c);
  exp2d(b,b);
  add21d(b,-2,c);
  test2(20,c);

  dto2d(20,a);
  chs2d(a, a);
  nroot2d(a,5,b);
  pow2di(b,5,c);
  sub2d(c,a,c);
  test2(21,c);

  dto2d(20,a);
  nroot2d(a,4,b);
  pow2di(b,4,c);
  sub2d(c,a,c);
  test2(22,c);

  for k:=1 to 99 do begin
    dto2d(k,a);
    {
    exp2d(a,b);
    ln2d(b,c);
    }
    ln2d(a,b);
    exp2d(b,c);
    sub2d(c,a,c);
    test2(k+100,c);
  end;

  for k:=78 to 99 do begin
    n := 1 + k mod 10;
    dto2d(k,a);
    dto2d(n,d);
    nroot2d(a,n,b);
    pow2d(b,d,c);
    sub2d(c,a,c);
    test2(k+200,c);
  end;
  inc(total_cnt, cnt);
  if fail>0 then begin
    writeln('*** failed');
    inc(total_failed,fail);
  end
  else writeln(cnt:4, ' tests OK')
end;


end.

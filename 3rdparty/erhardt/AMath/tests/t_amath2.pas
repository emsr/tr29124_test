{Part 2 of regression test for AMath unit,   (c) 2018  W.Ehrhardt}
{Tests for double-extended functions}

unit t_amath2;

{$i STD.INC}

{$ifdef BIT16}
{$N+}
{$endif}

interface

procedure test_xdivrem;
procedure test_fma_x;
procedure test_misc2x;


implementation

uses
  amath, t_amath0;


{---------------------------------------------------------------------------}
procedure test_xdivrem;
var
  i,k: integer;
  failed: longint;
  a,b,x,y,f,eps: extended;
begin
  writeln('Function xdivrem');
  failed := 0;
  {$ifdef FPC}
    eps := eps_d;        {FPC242 ... 264}
  {$else}
    eps := 0.5*eps_d;
  {$endif}
  for i:=-50 to 50 do begin
    if i<>0 then begin
      for k:=1 to 1000 do begin
        a := i;
        b := k;
        xdivrem(a,b,x,y);
        f := b*x + y;
        if abs(a-f) > eps*abs(a) then begin
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
procedure test_fma_x;
var
  x,y,f: double;
begin
  write('Test fma_x ... ');
  x := ldexp(1,-40);
  y := fma_x(1+x,1-x,-1);
  f := -x*x;
  if reldevx(y,f) > EPS then begin
    writeln('*** failed: ', reldevx(y,f));
    inc(total_failed);
  end
  else  writeln('OK.');
  inc(total_cnt);
end;


{---------------------------------------------------------------------------}
procedure test_misc2x;
var
  a,b,c,d: ext2;
  fail, cnt, k, n: integer;
const
  EPS = 5e-36;
const
  lnpih:  THexExtW = ($E85F,$3D0D,$8247,$9286,$3FFF);  {1.1447298858494001742}
  lnpil:  THexExtW = ($B668,$7BC0,$9395,$A06A,$BFBE);  {-3.3969475904465880268E-20}
  exppih: THexExtW = ($CCF6,$8499,$2375,$B920,$4003);  {2.3140692632779269005E+1}
  exppil: THexExtW = ($94D3,$AB76,$89F9,$DA74,$3FC2);  {7.4015511037705346451E-19}

  procedure test1(num: integer; x: ext2; i: longint);
  begin
    inc(cnt);
    if (x.h<>i) or (x.l<>0.0) then begin
      writeln('Test ', num, ' failed: ',i,' <> ',x.h:26,x.l:20);
    end;
  end;

  procedure test2(num: integer; x: ext2);
  begin
    inc(cnt);
    if abs(x.h+x.l) > EPS then begin
      inc(fail);
      writeln('Test ', num, ' failed: ',x.h:26,x.l:26);
    end;
  end;

begin
  writeln('Test other ext2 functions ');
  cnt := 0;
  fail := 0;

  xxto2x(1, 1e-70,a);
  ceil2x(a,b);
  test1(1, b, 2);
  floor2x(a,b);
  test1(2, b, 1);
  trunc2x(a,b);
  test1(3, b, 1);

  xxto2x(1, -1e-70,a);
  ceil2x(a,b);
  test1(4, b, 1);
  floor2x(a,b);
  test1(5, b, 0);
  trunc2x(a,b);
  test1(6, b, 0);

  xxto2x(-1, 1e-70,a);
  ceil2x(a,b);
  test1(7, b, 0);
  floor2x(a,b);
  test1(8, b, -1);
  trunc2x(a,b);
  test1(9, b, 0);

  xxto2x(-1, -1e-70,a);
  ceil2x(a,b);
  test1(10, b, -1);
  floor2x(a,b);
  test1(11, b, -2);
  trunc2x(a,b);
  test1(12, b, -1);

  xto2x(1,a);
  exp2x(a,b);
  exp1_2x(c);
  sub2x(c,b,c);
  test2(13,c);

  exp1_2x(a);
  ln2x(a,b);
  add21x(b,-1,c);
  test2(14, c);

  pi2x(a);
  ln2x(a,b);
  hhto2x(lnpih, lnpil, d);
  sub2x(b,d,c);
  test2(15,c);
  exp2x(d,b);
  sub2x(b,a,c);
  test2(16,c);

  pi2x(a);
  exp2x(a,b);
  hhto2x(exppih, exppil, d);
  sub2x(b,d,c);
  test2(17,c);
  ln2x(d,b);
  sub2x(b,a,c);
  test2(18,c);

  xto2x(2,a);
  ln2x(a,b);
  ln2_2x(c);
  sub2x(b,c,c);
  test2(19,c);
  exp2x(b,b);
  add21x(b,-2,c);
  test2(20,c);

  xto2x(20,a);
  chs2x(a, a);
  nroot2x(a,5,b);
  pow2xi(b,5,c);
  sub2x(c,a,c);
  test2(21,c);

  xto2x(20,a);
  nroot2x(a,4,b);
  pow2xi(b,4,c);
  sub2x(c,a,c);
  test2(22,c);

  for k:=1 to 99 do begin
    xto2x(k,a);
    {
    exp2x(a,b);
    ln2x(b,c);
    }
    ln2x(a,b);
    exp2x(b,c);
    sub2x(c,a,c);
    test2(k+100,c);
  end;

  for k:=78 to 99 do begin
    n := 1 + k mod 10;
    xto2x(k,a);
    xto2x(n,d);
    nroot2x(a,n,b);
    pow2x(b,d,c);
    sub2x(c,a,c);
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

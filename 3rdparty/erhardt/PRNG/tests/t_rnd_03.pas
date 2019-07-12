{Simple test for pasrand unit, we 12.2017}

program t_rnd_03;

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
  hrtimer, pasrand;

const
  LOOPS  = 10;

var
  ctx: pasrand_ctx;


{---------------------------------------------------------------------------}
function GenerateCycles: longint;
var
  i: integer;
  cyc0, cyc1, cyc2: comp;
  t1,t2,c1,c2: longint;
begin
  c1 := MaxLongint;
  c2 := MaxLongint;
  for i:=1 to LOOPS do begin
    ReadTSC(cyc0);
    pasrand_next(ctx);
    ReadTSC(cyc1);
    pasrand_next(ctx);
    pasrand_next(ctx);
    pasrand_next(ctx);
    pasrand_next(ctx);
    pasrand_next(ctx);
    ReadTSC(cyc2);
    t2 := round(cyc2-cyc1);
    t1 := round(cyc1-cyc0);
    if t1<c1 then c1 := t1;
    if t2<c2 then c2 := t2;
  end;
  GenerateCycles := (c2-c1+1) shr 2;
end;


var
  CBlk: longint;
  CPB : double;
begin
  writeln('Test for pasrand unit  (c) 2017 W.Ehrhardt');
  writeln(' pasrand selftest: ',pasrand_selftest);
  fillchar(ctx, sizeof(ctx),0);
  pasrand_init(ctx, 0);
  CBlK := GenerateCycles;
  CPB  := CBlk/4.0;
  writeln('   CPU Frequency: ', CPUFrequency/1E6:1:1);
  writeln('        Generate: ', CBlk);
  writeln('     Cycles/Byte: ', CPB:1:1);
  writeln('            MB/s: ', CPUFrequency/CPB/1E6:1:3);
end.

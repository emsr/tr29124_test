{Simple test for aesr unit, we May 2005}

program t_rnd_73;

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
  hrtimer, aesr;

const
  LOOPS  = 10;

var
  ctx: aesr_ctx;



{---------------------------------------------------------------------------}
procedure aesr_generate(var ctx: aesr_ctx);
  {-generate next 256 result values}
var
  i: integer;
begin
  for i:=0 to 3 do aesr_next(ctx);
end;


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
    aesr_generate(ctx);
    ReadTSC(cyc1);
    aesr_generate(ctx);
    aesr_generate(ctx);
    aesr_generate(ctx);
    aesr_generate(ctx);
    aesr_generate(ctx);
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
  writeln('Test for aesr unit  (c) 2005 W.Ehrhardt');
  writeln('   aesr selftest: ',aesr_selftest);
  fillchar(ctx, sizeof(ctx),0);
  aesr_init(ctx, 0);
  CBlK := GenerateCycles;
  CPB  := CBlk/16;
  writeln('   CPU Frequency: ', CPUFrequency/1E6:1:1);
  writeln('        Generate: ', CBlk);
  writeln('     Cycles/Byte: ', CPB:1:1);
  writeln('            MB/s: ', CPUFrequency/CPB/1E6:1:3);
end.

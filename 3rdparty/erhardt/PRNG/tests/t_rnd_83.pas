{Simple test for salsar unit, we Apr.2006}

program t_rnd_83;

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
  hrtimer, salsar;

const
  LOOPS  = 10;

var
  ctx: salsar_ctx;



{---------------------------------------------------------------------------}
procedure salsar_generate(var ctx: salsar_ctx);
  {-generate next SR_KSB_SIZE result values}
var
  i: integer;
begin
  for i:=1 to SR_KSB_SIZE do salsar_next(ctx);
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
    salsar_generate(ctx);
    ReadTSC(cyc1);
    salsar_generate(ctx);
    salsar_generate(ctx);
    salsar_generate(ctx);
    salsar_generate(ctx);
    salsar_generate(ctx);
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
  writeln('Speed test for salsar unit  (c) 2006 W.Ehrhardt');
  writeln(' salsar selftest: ', salsar_selftest);
  writeln('   No. of rounds: ', salsar_get_rounds);
  fillchar(ctx, sizeof(ctx),0);
  salsar_init(ctx, 0);
  CBlK := GenerateCycles;
  CPB  := CBlk/(4*SR_KSB_SIZE);
  writeln('   CPU Frequency: ', CPUFrequency/1E6:1:1);
  writeln('        Generate: ', CBlk);
  writeln('     Cycles/Byte: ', CPB:1:1);
  writeln('            MB/s: ', CPUFrequency/CPB/1E6:1:3);
end.

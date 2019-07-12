{Simple test for salsar unit, we Apr.2006}

program t_rnd_80;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  salsar;

var
  stream: array[0..127] of longint;
  sarr  : salsar_sarr absolute stream;
  ctx   : salsar_ctx;
  i     : integer;
const
  rounds: array[0..2] of word = (8,12,20);
  {Test data from Set 1, vector #9; values for 128 bit key size}
  test  : array[0..2] of longint = (longint($e7a698c4),longint($55f454e3),longint($0c141552));
(*
  {Test values for 256 bit key size}
  test  : array[0..2] of longint = (longint($7026acd4),longint($7caf2580),longint($9e9b60b0));
*)
begin
  writeln('Simple test for salsar unit     (c) 2006 W.Ehrhardt');
  for i:=0 to 2 do begin
    salsar_set_rounds(rounds[i]);
    writeln('Number of rounds: ', rounds[i]);
    writeln('salsar self test: ', salsar_selftest);
    fillchar(stream,sizeof(stream),0);
    sarr[0] := $4000;
    salsar_inita(ctx, sarr);
    salsar_read(ctx, @stream, sizeof(stream));
    writeln('salsar_read test: ', stream[127]=test[i]);
  end;
end.

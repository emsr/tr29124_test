program t_ldemo;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

{$i std.inc}

uses
{$IFnDEF FPC}
{$ELSE}
  Interfaces,
{$ENDIF}
  Forms,
  t_ldemou in 't_ldemou.pas' {Form1};

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.

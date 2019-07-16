program t_ddemo;

{$i std.inc}

uses
  {$ifdef UNIT_SCOPE}
    vcl.Forms,
  {$else}
    Forms,
  {$endif}
  t_ddemou in 't_ddemou.pas' {Form1};

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.

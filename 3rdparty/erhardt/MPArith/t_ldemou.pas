{***************************************************************
 *
 * Unit Name: t_ddemou
 * Purpose  : Delphi MPArith demo program (calculate Mersenne primes)
 * Author   : W.Ehrhardt Sep.2005
 * History  : V 0.1  Sep.2005
                0.2  Dec.2012  XE3
 *
 ****************************************************************}

unit t_ldemou;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

interface

{$i std.inc}

uses
  {$ifdef UNIT_SCOPE}
    winapi.Windows, winapi.Messages, system.SysUtils,
    system.Classes, vcl.Graphics, vcl.Controls, vcl.Forms, vcl.Dialogs,
    vcl.StdCtrls, vcl.Buttons, vcl.ExtCtrls,
  {$else}
{$IFnDEF FPC}
  Windows,
{$ELSE}
  LCLIntf, LCLType, LMessages,
{$ENDIF}
  Messages, SysUtils, Classes, Graphics, Controls, Forms,
    Dialogs, StdCtrls, Buttons, ExtCtrls,
  {$endif}
  mp_types, mp_numth;

type
  TForm1 = class(TForm)
    Panel1: TPanel;
    Memo1: TMemo;
    BitBtn_Start: TSpeedButton;
    BitBtn_Stop: TSpeedButton;
    Label_Current: TLabel;
    procedure FormShow(Sender: TObject);
    procedure BitBtn_StartClick(Sender: TObject);
    procedure BitBtn_StopClick(Sender: TObject);
    procedure FormCloseQuery(Sender: TObject; var CanClose: Boolean);
  private
    { Private declarations }
    s: integer;              //Current index for Mersenne(s)
    abort: boolean;          //Stop button pressed
    running: boolean;        //Calculation is already running
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation

{$IFnDEF FPC}
  {$R *.dfm}
{$ELSE}
  {$R *.lfm}
{$ENDIF}


{---------------------------------------------------------------------------}
procedure TForm1.FormShow(Sender: TObject);
  {-Show mpint version and program purpose}
begin
  Memo1.Lines.Add('Test of MPArith '+ MP_VERSION + '   (c) W.Ehrhardt 2005-2018');
  Memo1.Lines.Add('Calculate Mersenne primes using Lucas-Lehmer test');
  Memo1.Lines.Add('');
  s:=1;
  running := false;
end;


{---------------------------------------------------------------------------}
procedure TForm1.BitBtn_StartClick(Sender: TObject);
  {-Loop with mp_isMersennePrime, update memo and current index}
begin
  if running then exit;
  {Set control flags}
  running := true;
  abort := false;
  BitBtn_Start.Enabled := false;
  while (s<MaxMersenne) and not abort do begin
    Label_Current.Caption := 'Current: '+IntToStr(s);
    {Give a chance to update the label and check for stop button}
    Application.ProcessMessages;
    if abort then break;
    if mp_isMersennePrime(s) then begin
      Memo1.Lines.Add('Prime: 2^'+IntToStr(s)+'-1');
    end;
    inc(s);
  end;
  BitBtn_Start.Enabled := true;
  running := false;
end;


{---------------------------------------------------------------------------}
procedure TForm1.BitBtn_StopClick(Sender: TObject);
  {-signal abort if stop button clicked}
begin
  abort:=true;
end;


{---------------------------------------------------------------------------}
procedure TForm1.FormCloseQuery(Sender: TObject; var CanClose: Boolean);
  {-signal abort if application wants to close}
begin
  if running and not abort then CanClose := false;
  abort := true;
end;


end.

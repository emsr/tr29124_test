
{---------------------------------------------------------------------------}
{------------ mp include file for taus88 random number generator -----------}
{---------------------------------------------------------------------------}


uses taus88;

var
  mp_ctx: taus88_ctx;      {Global taus88 context}


{---------------------------------------------------------------------------}
function  mp_random_byte: byte;
  {-Returns a random byte}
begin
  mp_random_byte := taus88_word(mp_ctx) and $FF;
end;


{---------------------------------------------------------------------------}
function mp_random_long: longint;
  {-Returns a random positive longint}
begin
  mp_random_long := taus88_long(mp_ctx)
end;


{---------------------------------------------------------------------------}
function mp_random_int: longint;
  {-Returns a random signed longint}
begin
  taus88_next(mp_ctx);
  mp_random_int := mp_ctx.nr;
end;


{---------------------------------------------------------------------------}
function mp_random_digit: mp_digit;
  {-Returns a random mp_digit}
begin
  taus88_next(mp_ctx);
  mp_random_digit := mp_digit(mp_ctx.nr and MP_DIGIT_MAX);
end;


{---------------------------------------------------------------------------}
procedure mp_random_read(dest: pointer; len: word);
  {-Read len bytes from the PRNG to dest}
begin
  taus88_read(mp_ctx, dest, len);
end;


{---------------------------------------------------------------------------}
procedure mp_random_seed(const seed: array of longint);
  {-Initialize PRNG with array of longint}
var
  s1,s2,s3: longint;
begin
  s1 := seed[0];
  if high(seed)>0 then s2 := seed[1] else s2 := 0;
  if high(seed)>1 then s3 := seed[2] else s3 := 0;
  taus88_init3(mp_ctx,s1,s2,s3);
end;


{---------------------------------------------------------------------------}
procedure mp_random_randomize;
  {-Initialize PRNG via randomize/TSC}
{$ifdef MPC_UseTSC}
var
  seed: array[0..2] of longint;
{$endif}
begin
  randomize;
  {$ifdef MPC_UseTSC}
    seed[0] := RandSeed;
    {$ifdef BASM16}
      asm
        db  $0f,$31
        db  $66;  mov word ptr [seed+4], ax
        db  $66;  mov word ptr [seed+8], dx
      end;
    {$else}
      asm
         dw  $310F
         mov dword ptr [seed+4], eax
         mov dword ptr [seed+8], edx
      end;
    {$endif}
    mp_random_seed(seed);
  {$else}
    taus88_init0(mp_ctx);
  {$endif}
end;


{---------------------------------------------------------------------------}
function  mp_random_word: word;
  {-Returns a random word}
begin
  mp_random_word := taus88_word(mp_ctx);
end;


{---------------------------------------------------------------------------}
{------------ mp include file for ISAAC random number generator ------------}
{---------------------------------------------------------------------------}


uses isaac;

var
  mp_ctx: isaac_ctx;      {Global ISAAC context}


{---------------------------------------------------------------------------}
function  mp_random_byte: byte;
  {-Returns a random byte}
begin
  mp_random_byte := isaac_word(mp_ctx) and $FF;
end;


{---------------------------------------------------------------------------}
function mp_random_long: longint;
  {-Returns a random positive longint}
begin
  mp_random_long := isaac_long(mp_ctx)
end;


{---------------------------------------------------------------------------}
function mp_random_int: longint;
  {-Returns a random signed longint}
begin
  isaac_next(mp_ctx);
  mp_random_int := mp_ctx.nextres;
end;


{---------------------------------------------------------------------------}
function mp_random_digit: mp_digit;
  {-Returns a random mp_digit}
begin
  isaac_next(mp_ctx);
  mp_random_digit := mp_digit(mp_ctx.nextres and MP_DIGIT_MAX);
end;



{---------------------------------------------------------------------------}
procedure mp_random_read(dest: pointer; len: word);
  {-Read len bytes from the PRNG to dest}
begin
  isaac_read(mp_ctx, dest, len);
end;


{---------------------------------------------------------------------------}
procedure mp_random_seed(const seed: array of longint);
  {-Initialize PRNG with array of longint}
begin
  isaac_inita(mp_ctx, seed, succ(High(seed)));
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
    isaac_init0(mp_ctx);
  {$endif}
end;


{---------------------------------------------------------------------------}
function  mp_random_word: word;
  {-Returns a random word}
begin
  mp_random_word := isaac_word(mp_ctx);
end;

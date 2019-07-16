{Basic 16 bit routines}

{---------------------------------------------------------------------------}
function bitsize32(a: longint): integer; assembler;
  {-Return the number of bits in a (index of highest bit), 0 if no bit is set}
asm
   db $66;  mov ax,word ptr [a]
   db $66;  or  ax,ax
            je  @x
            db  $66,$0F,$BD,$C0   {bsr  eax,eax }
            inc ax;
   @x:
end;


{---------------------------------------------------------------------------}
function __gcd32: longint; near; assembler;
  {-Calculate GCD of unsigned (A,B) in (eax,edx)}
asm
{done if A=B, otherwise if A<B then swap A and B}
     db $66;  cmp   ax,dx
              jz    @@x
              jae   @@1
     db $66;  xchg  ax,dx

{here eax >= edx. Calculate odd parts a,b with A=a*2^e(a), B=b*2^e(b)}
@@1: db $66,$0F,$BC,$CA       {bsf ecx,edx; if B=0 return A}
              jz    @@x

     db $66,$0F,$BC,$D8       {bsf ebx,eax; ebx=e(a), A cannot be zero}
     db $66;  shr   dx, cl    {edx=b}
     db $66;  xchg  bx, cx
     db $66;  shr   ax, cl    {eax=a}

     db $66;  cmp   bx,cx     {ebx = e = min(e(a),e(b)}
              jb    @@2
     db $66;  mov   bx,cx

@@2: db $66;  cmp   ax,dx     {compare a and b}
              jz    @@4       {done if equal}

{in the main loop both a and b are always odd}
{therefore for |a-b| is even and non-zero}

@@3: db $66;  mov   si, ax    {eax=a, edx=b}
     {calculate max(a,b) and min(a,b) without branches}
     {see H.S.Warren, Hacker's Delight, Revision 1/4/07}
     {http://www.hackersdelight.org/revisions.pdf}
     db $66;  sub   si, dx    {esi=a-b}
     db $66;  sbb   cx, cx    {if a>=b then ecx=0 else ecx=-1}
     db $66;  and   si, cx    {if a>=b then esi=0 else esi=a-b}
     db $66;  add   dx, si    {if a>=b then edx=b else edx=a, i.e. edx=min(a,b)}
     db $66;  sub   ax, si    {if a>=b then eax=a else eax=a-(a-b)=b=max(a,b)}
     db $66;  sub   ax, dx    {a'=max(a,b)-min(a,b), b'=min(a,b)}
     db $66,$0F,$BC,$C8       {bsf ecx,eax: a'=|a-b| is even, divide out powers of 2}
     db $66;  shr   ax, cl
     db $66;  cmp   ax, dx    {compare new a and new b}
              jnz   @@3       {and repeat loop if not equal}

@@4: db $66;  mov   cx, bx    {shift by initial common exponent e}
     db $66;  shl   ax, cl
@@x: db $66;  mov   dx,ax     {return GCD in (dx:ax)}
     db $66;  shr   dx,16
end;


{---------------------------------------------------------------------------}
function gcd32u(A, B: longint): longint; assembler;
  {-Calculate GCD of two longints (DWORD interpretation)}
asm
     db $66;  mov   ax,word ptr [A]
     db $66;  mov   dx,word ptr [B]
              call  __gcd32
end;


{---------------------------------------------------------------------------}
function gcd32(A, B: longint): longint; assembler;
  {-Calculate GCD of two longints}
asm
     db $66;  mov   ax,word ptr [A]
     db $66;  mov   dx,word ptr [B]
     db $66;  and   ax,ax
              jns   @@1
     db $66;  neg   ax
@@1: db $66;  and   dx,dx
              jns   @@2
     db $66;  neg   dx
@@2:          call  __gcd32
end;


{---------------------------------------------------------------------------}
function mulmod32(a,b,n: longint): longint; assembler;
  {-Return a*b mod n, assumes n>0, a,b>=0}
asm
  db $66; mov   ax,word ptr [a]
  db $66; mov   dx,word ptr [b]
  db $66; mov   cx,word ptr [n]
  db $66; mul   dx
  db $66; div   cx
  db $66; mov   ax,dx   {return result in (dx:ax)}
  db $66; shr   dx,16
end;


{---------------------------------------------------------------------------}
function invmod32(a,b: longint): longint;
  {-Return a^-1 mod b, b>1. Result is 0 if gcd(a,b)<>1 or b<2}
var
  g,i: longint;
begin
  if (b>1) and (a<>0) then begin
    {Use extended GCD to calculate u1*a + u2*b = u3 = gcd(a.b)  }
    {If u3 = 1, then u1 = a^-1 mod b. u2 will not be calculated.}
    {Notation from Knuth [3] Algorithm X. u3 and v3 will be >=0 }
    {and |u1| <= b, |v1| <= b, see e.g. Shoup [29], Theorem 4.3.}
    {u1 = ebx = 1  }
    {u3 = ecx = |a|}
    {v1 = esi = 0  }
    {v3 = edi = b  }
    asm
           db $66;   mov    di,word ptr [b]
           db $66;   mov    ax,word ptr [a]
           db $66;   cwd    {cdq!}
           db $66;   xor    ax,dx
           db $66;   sub    ax,dx
           db $66;   mov    cx,ax   {ecx=u3=abs(a)}
           db $66;   sub    si,si
           db $66;   sub    bx,bx
                     inc    bx
      @@1: db $66;   or     di,di   {done if v3=0}
                     jz     @@2
           db $66;   mov    ax,cx
           db $66;   sub    dx,dx
           db $66;   div    di      {eax=q=u3 div v3, edx=u3 mod v3}
           db $66;   mov    cx,di   {u3' = v3}
           db $66;   mov    di,dx   {v3' = u3 mod v3}
           db $66;   imul   si
           db $66;   sub    bx,ax   {ebx=u1-q*v1}
           db $66;   xchg   bx,si
                     jmp    @@1
      @@2: db $66;   mov    word ptr [g], cx
           db $66;   mov    word ptr [i], bx
    end;
    if g=1 then begin
      {gcd(a,b)=1, so inverse exists: do some sign related adjustments.}
      if i<0 then inc(i,b);
      if (a<0) and (i<>0) then invmod32 := b-i else invmod32 := i;
    end
    else invmod32 := 0;
  end
  else invmod32 := 0;
end;


{---------------------------------------------------------------------------}
function add32_ovr(x,y: longint; var z: longint): boolean;
  {-Add z=x+y with overflow detection}
var
  overflow: boolean;
begin
  asm
            sub  dx,dx
    db $66; mov  ax,word ptr [x]
    db $66; add  ax,word ptr [y]
            jno  @@1
            inc  dx
    @@1:    mov  [overflow],dl
            les  di, [z]
    db $66; mov  [es:di],ax
  end;
  add32_ovr := overflow;
end;


{---------------------------------------------------------------------------}
function fLeftShiftAdd(w: longint; d: word): longint;
  {-Return (w shl DIGIT_BIT) + d}
inline(
  $66/$2B/$D2/             {sub    edx,edx      }
  $5A/                     {pop    dx           }
  $66/$58/                 {pop    eax          }
  $66/$C1/$E0/<DIGIT_BIT/  {shl    eax,DIGIT_BIT}
  $66/$03/$C2/             {add    eax,edx      }
  $66/$8B/$D0/             {mov    edx,eax      }
  $66/$C1/$EA/$10);        {shr    edx,16       }


{---------------------------------------------------------------------------}
procedure LeftShiftAdd(var w: longint; d: word);
  {-Calculate w := (w shl DIGIT_BIT) + d}
inline(
  $66/$2B/$D2/             {sub    edx,edx      }
  $5A/                     {pop    dx           }
  $5E/                     {pop    si           }
  $07/                     {pop    es           }
  $66/$26/$8B/$04/         {mov    eax,es:[si]  }
  $66/$C1/$E0/<DIGIT_BIT/  {shl    eax,DIGIT_BIT}
  $66/$03/$C2/             {add    eax,edx      }
  $66/$26/$89/$04);        {mov    es:[si],eax  }


{--------------------------------------------------------------------------}
function HeapFunc(Size: word): integer; far;
  {-Forces nil return values instead of runtime error if out of memory}
begin
  if size>0 then HeapFunc := 1;
end;


{---------------------------------------------------------------------------}
function  IAlloc(lsize: longint): pointer;
  {-Allocate heap, return nil if error, no diagnostics}
var
  p, SaveHeapError : pointer;
type
  LH = packed record L,H: word; end;
begin
  if LH(lsize).H<>0 then IAlloc := nil
  else begin
    SaveHeapError := HeapError;
    HeapError := @HeapFunc;
    getmem(p, lsize);
    HeapError := SaveHeapError;
    IAlloc  := p;
  end;
end;


{---------------------------------------------------------------------------}
function mp_realloc(p: pointer; oldsize, newsize: longint): pointer;
  {-Reallocate heap to new size, normally new pointer will be <> p.}
  { If newsize>oldsize the new allocated space is zerofilled}
var
  tmp: pointer;
begin
  if oldsize=newsize then begin
    mp_realloc := p;
    exit;
  end;
  tmp := mp_getmem(newsize);
  mp_realloc := tmp;
  if tmp=nil then exit;
  if newsize>oldsize then begin
    {copy and zero fill}
    move(p^,tmp^,oldsize);
    inc(Ptr2Inc(tmp),oldsize);
    fillchar(tmp^,newsize-oldsize,0);
  end
  else begin
    {newsize <= oldsize: copy only}
    move(p^,tmp^,newsize);
  end;
  if mp_clearzero then fillchar(p^, oldsize, 0);
  mp_freemem(p,oldsize);
end;

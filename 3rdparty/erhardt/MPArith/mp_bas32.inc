{Basic 32 bit routines}

{---------------------------------------------------------------------------}
function bitsize32(a: longint): integer; assembler;
  {-Return the number of bits in a (index of highest bit), 0 if no bit is set}
asm
  {$ifdef LoadArgs}
    mov eax,[a]
  {$endif}
    or    eax, eax
    je    @x
    bsr   eax, eax
    inc   eax
  @x:
end;


{---------------------------------------------------------------------------}
function __gcd32: longint; assembler;
  {-Calculate GCD of unsigned (A,B) in (eax,edx)}
asm
        push  ebx
        {done if A=B, otherwise if A<B then swap A and B}
        cmp   eax,edx
        jz    @@x
        jae   @@1
        xchg  eax,edx

{here eax >= edx. Calculate odd parts a,b with A=a*2^e(a), B=b*2^e(b)}
@@1:    bsf   ecx, edx         {if B=0 return A}
        jz    @@x

        bsf   ebx, eax         {ebx=e(a), A cannot be zero}
        shr   edx, cl          {edx=b}
        xchg  ebx, ecx
        shr   eax, cl          {eax=a}

        cmp   ebx,ecx          {ebx = e = min(e(a),e(b)}
        jb    @@2
        mov   ebx,ecx

@@2:    cmp   eax,edx          {compare a and b}
        jz    @@4              {done if equal}

{in the main loop both a and b are always odd}
{therefore for |a-b| is even and non-zero}
        push  esi
@@3:    mov   esi, eax         {eax=a, edx=b}
        {calculate max(a,b) and min(a,b) without branches}
        {see H.S.Warren, Hacker's Delight, Revision 1/4/07}
        {http://www.hackersdelight.org/revisions.pdf}
        sub   esi, edx         {esi=a-b}
        sbb   ecx, ecx         {if a>=b then ecx=0 else ecx=-1}
        and   esi, ecx         {if a>=b then esi=0 else esi=a-b}
        add   edx, esi         {if a>=b then edx=b else edx=a, i.e. edx=min(a,b)}
        sub   eax, esi         {if a>=b then eax=a else eax=a-(a-b)=b=max(a,b)}
        sub   eax, edx         {a'=max(a,b)-min(a,b), b'=min(a,b)}
        bsf   ecx, eax         {a'=|a-b| is even, divide out powers of 2}
        shr   eax, cl
        cmp   eax, edx         {compare new a and new b}
        jnz   @@3              {and repeat loop if not equal}

        pop   esi

@@4:    mov   ecx, ebx         {shift by initial common exponent e}
        shl   eax, cl
@@x:    pop   ebx

end;


{---------------------------------------------------------------------------}
function  gcd32u(A, B: longint): longint; assembler;
  {-Calculate GCD of two longints (DWORD interpretation)}
asm
  {$ifdef LoadArgs}
        mov   eax,[A]
        mov   edx,[B]
  {$endif}
        call  __gcd32
end;


{---------------------------------------------------------------------------}
function gcd32(A, B: longint): longint; assembler;
  {-Calculate GCD of two longints}
asm
  {$ifdef LoadArgs}
       mov   eax,[A]
       mov   edx,[B]
  {$endif}
       and   eax,eax
       jns   @@1
       neg   eax
@@1:   and   edx,edx
       jns   @@2
       neg   edx
@@2:   call  __gcd32
end;


{---------------------------------------------------------------------------}
function mulmod32(a,b,n: longint): longint;
  {-Return a*b mod n, assumes n>0, a,b>=0}
begin
  asm
    mov   eax,[a]
    mov   edx,[b]
    mov   ecx,[n]
    mul   edx
    div   ecx
  {$ifdef FPC}
    mov   [a],edx
  {$else}
    mov   [@result],edx
  {$endif}
  end;
  {$ifdef FPC}
    mulmod32 := a;
  {$endif}
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
    asm    push   ebx
           push   esi
           push   edi
           mov    edi,[b]
           mov    eax,[a]
           cdq
           xor    eax,edx
           sub    eax,edx
           mov    ecx,eax     {ecx=u3=abs(a)}
           sub    esi,esi
           mov    ebx,1
      @@1: or     edi,edi     {done if v3=0}
           jz     @@2
           mov    eax,ecx
           sub    edx,edx
           div    edi         {eax=q=u3 div v3, edx=u3 mod v3}
           mov    ecx,edi     {u3' = v3}
           mov    edi,edx     {v3' = u3 mod v3}
           imul   esi
           sub    ebx,eax     {ebx=u1-q*v1}
           xchg   ebx,esi
           jmp    @@1
      @@2: mov    [g], ecx
           mov    [i], ebx
           pop    edi
           pop    esi
           pop    ebx
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
            sub  edx,edx
            mov  eax,[x]
            add  eax,[y]
            jno  @@1
            inc  edx
    @@1:    mov  [overflow],dl
            mov  edx,[z]
            mov  [edx],eax
  end;
  add32_ovr := overflow;
end;


{$ifdef FPC}

{---------------------------------------------------------------------------}
function IAlloc(lsize: longint): pointer;
  {-Allocate heap > 64K, return nil if error, no diagnostics}
var
  p: pointer;
  sh: boolean;
begin
  sh := ReturnNilIfGrowHeapFails;
  ReturnNilIfGrowHeapFails := true;
  getmem(p, lsize);
  ReturnNilIfGrowHeapFails := sh;
  IAlloc:= p;
end;

{$else}

{---------------------------------------------------------------------------}
function IAlloc(lsize: longint): pointer;
  {-Allocate heap > 64K, return nil if error, no diagnostics}
var
  p: pointer;
begin
  try
    getmem(p, lsize);
  except
    p := nil;
  end;
  IAlloc := p;
end;

{$endif FPC}

{---------------------------------------------------------------------------}
function mp_realloc(p: pointer; oldsize, newsize: longint): pointer;
  {-Reallocate heap to new size, if newsize>oldsize the new allocated space is zerofilled}
var
  tmp: pointer;
begin
  if oldsize=newsize then begin
    mp_realloc := p;
    exit;
  end;
  tmp := p;
  ReallocMem(tmp, newsize);
  mp_realloc := tmp;
  if tmp<>nil then begin
    {$ifdef MPC_Diagnostic}
      inc(mp_memstat.MemDiff, newsize);
      dec(mp_memstat.MemDiff, oldsize);
    {$endif}
    if newsize>oldsize then begin
      {zero fill new part}
      inc(Ptr2Inc(tmp),oldsize);
      fillchar(tmp^,newsize-oldsize,0);
    end;
  end;
end;

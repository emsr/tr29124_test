
;          Copyright Christopher Kormanyos 2002 - 2009.
; Distributed under the Boost Software License, Version 1.0.
;    (See accompanying file LICENSE_1_0.txt or copy at
;          http://www.boost.org/LICENSE_1_0.txt)

; This is a legacy file which is not used
.686
.XMM


.MODEL FLAT, C

.CODE


__scale_32_SSE  DQ 2882303760
__ef_elem_mask  DQ 100000000


; void mul_loop_uv(const UINT32* const u, const UINT32* const v, UINT32* const w, const INT32 jmax,  const INT32 imax)
; {
;   for(INT32 j = jmax - 1; j >= 0; j--)
;   {
;     UINT32 carry = static_cast<UINT32>(0u);
;
;     INT32 i = (imax - 1) - j;
;
;     if(u[j] != static_cast<UINT32>(0u))
;     {
;       for( ; i >= 0; i--)
;       {
;         // The following lines are equivalent to:
;         //   carry          = t / mask
;         //   w[(i + j) + 1] = t % mask
;
;         const UINT46 t = static_cast<UINT64>(static_cast<UINT64>(u[j] * static_cast<UINT64>(v[i])) + static_cast<UINT32>(w[(i + j) + 1] + carry));
;         carry          = static_cast<UINT32>(t / static_cast<UINT32>(efx::e_float::ef_elem_mask));
;         w[(i + j) + 1] = static_cast<UINT32>(t - static_cast<UINT64>(carry * static_cast<UINT64>(efx::e_float::ef_elem_mask)));
;       }
;     }
; 
;     w[(i + j) + 1] = carry;
;   }
; }


mul_loop_uv PROC

  pushad

  mov ebx, [esp + 030h]                       ; ebx : the value of jmax is in ebx
  dec ebx                                     ; ebx : the value of jmax - 1 is in ebx
  mov ebp, [esp + 024h]                       ; ebp : the address &u[0] is in ebp

  loop_mul_outer:

    mov ecx, [esp + 034h]                     ; ecx: the value of imax is in ecx
    sub ecx, ebx                              ; ecx: the value of imax - j is in ecx
    dec ecx                                   ; ecx: the value of the index i = imax - 1 - j is in ecx

    xor edi, edi                              ; edi : edi is used for the carry;

    cmp dword ptr [ebp + 04h*ebx], 0h         ; if(u[j] != static_cast<UINT32>(0u))
    je loop_mul_outer_increment

    mov esi, ebx
    add esi, ecx                              ; eax : the value of i + j is in eax
    shl esi, 02h
    add esi, [esp + 02Ch]                     ; esi : esi points to the address of &w[(i + j)]

    loop_mul_inner:

      mov edx, [esp + 028h]                   ; edx : the address &v[0] is in edx
      mov eax, dword ptr [edx + 04h*ecx]      ; eax : the value of v[i] is in eax
      add edi, dword ptr [esi + 04h]          ; edi : w[(i + j) + 1] + carry

      mul dword ptr [ebp + 04h*ebx]           ; multiply v[i] with u[j]
      add eax, edi                            ; add (w[(i + j) + 1] + carry) to the product
      adc edx, 0h                             ; edx:eax : t = u[j] * v[i] + (w[(i + j) + 1] + carry)

      movd mm1, eax
      movd mm2, edx
      movd mm0, edx
      psllq mm0, 020h
      por mm0, mm1                            ; mm0 : t = u[j] * v[i] + (w[(i + j) + 1] + carry)

      ; Begin the division by starting with the 64*32=96 bit multiplication with the
      ; scale which is stored in mm4. Then shift the upper 64 bits of the 96-bit result
      ; to the right by 58 places (note: 58-32=26) to obtain the result of the division.

      pmuludq mm1, qword ptr [__scale_32_SSE] ; mm1  *= 2882303760;
      pmuludq mm2, qword ptr [__scale_32_SSE] ; mm2  *= 2882303760;
      psrlq mm1, 020h                         ; mm1 >>= 32u;
      paddq mm2, mm1                          ; mm2  += mm1;
      psrlq mm2, 01ah                         ; mm2 >>= 26u;

      ; Check if the division is complete by multiplying the product and subtracting
      ; this from the dividend. Correct the result if necessary if the result of
      ; the remainder is larger than or equal to 100000000. Note that the correction
      ; to the division can be no larger than 1 in the quotient.

      movd edi, mm2                           ; edi : the value of the quotient (the carry) is in edi
      pmuludq mm2, qword ptr [__ef_elem_mask] ; quotient *= 100000000;
      psubq mm0, mm2                          ; mm0 : the value of the remainder is in mm0
      movd dword ptr [esi + 04h], mm0         ; remainder : w[(i + j) + 1] = t % mask

      cmp dword ptr [esi + 04h], 100000000

      jl loop_mul_inner_increment

        inc edi
        sub dword ptr [esi + 04h], 100000000

      loop_mul_inner_increment:

      sub esi, 04h

      dec ecx
      cmp ecx, 0h
      jge loop_mul_inner

      mov dword ptr [esi + 04h], edi

    loop_mul_outer_increment:

    sub ebx, 1h
    test ebx, ebx
    jge loop_mul_outer

  popad

  emms

  ret
  mul_loop_uv ENDP


; UINT32 mul_loop_n(UINT32* const u, UINT32 n, const INT32 jm)
; {
;   UINT32 carry = static_cast<UINT32>(0u);
;
;   for(INT32 j = jmax - 1; j >= 1; j--)
;   {
;     const UINT64 t = static_cast<UINT64>(carry + static_cast<UINT64>(u[j] * static_cast<UINT64>(n)));
;     carry          = static_cast<UINT32>(t / static_cast<UINT32>(mask));
;     u[j]           = static_cast<UINT32>(t - static_cast<UINT64>(static_cast<UINT32>(mask) * carry));
;   }
;
;   return carry;
; }

mul_loop_n PROC

  pushad

  mov esi, 100000000

  mov edi, [esp + 2Ch]                  ; edi : the value of the jmax is in edi
  dec edi                               ; edi : the value of jmax - 1 is in edi

  mov ebx, [esp + 24h]                  ; the address &u[0] is in ebx
  movd mm7, dword ptr [esp + 28h]       ; the value of n is in mm7

  pxor mm0, mm0                         ; the value of the carry is in mm0

  loop_mul_n:

    movd mm1, dword ptr [ebx + 4h*edi]
    pmuludq mm1, mm7
    paddq mm0, mm1                      ; the value of t is in mm0
    movq mm3, mm0                       ; the value of t is in mm3

    movd eax, mm0                       ; divide the carry by the mask (100000000)
    psrlq mm0, 020h
    movd edx, mm0
    div esi                             ; divide by the mask (100000000)
    movd mm0, eax                       ; the value of the carry is, once again, in mm0
    mov dword ptr [ebx + 4h*edi], edx   ; store the result of (t % mask) in u[j - 1]

    sub edi, 1h
    test edi, edi
    jge loop_mul_n

  popad

  movd eax, mm0

  emms

  ret
  mul_loop_n ENDP


; UINT32 div_loop_n(UINT32* const u, UINT32 n, const INT32 jm)
; {
;   UINT32 prev = static_cast<UINT32>(0u);
;
;   for(INT32 j = static_cast<INT32>(0); j < jm; j++)
;   {
;     const UINT64 t = static_cast<UINT64>(u[j] + static_cast<UINT64>(prev * static_cast<UINT64>(mask)));
; 
;     u[j] = static_cast<UINT32>(t / n);
;     prev = static_cast<UINT32>(t - static_cast<UINT64>(static_cast<UINT32>(n) * static_cast<UINT64>(u[j])));
;   }
; 
;   return prev;
; }


div_loop_n PROC

  pushad

  movd mm5, dword ptr [__ef_elem_mask]  ; the value of the mask (100000000) is in mm5

  mov edi, [esp + 2Ch]                  ; the value of the jm is in edi

  mov ebx, [esp + 24h]                  ; the address &u[0] is in ebx
  movd mm7, dword ptr [esp + 28h]       ; the value of n is in mm7

  xor ecx, ecx
  pxor mm0, mm0                         ; the value of prev is in mm0

  loop_div_n:
  
    movd mm1, dword ptr [ebx + 4h*ecx]  ; the value of u[j] is in mm1
    pmuludq mm0, mm5
    paddq mm0, mm1                      ; the value of t is in mm0

    movd eax, mm0                       ; divide t by n
    psrlq mm0, 020h
    movd edx, mm0
    div dword ptr [esp + 28h]           ; n is in stack memory at [esp + 28h]
    mov dword ptr [ebx + 4h*ecx], eax   ; set u[j]
    movd mm0, edx                       ; the value of prev is in mm0

    inc ecx
    cmp ecx, edi
    jl loop_div_n

  popad

  movd eax, mm0

  emms

  ret
  div_loop_n ENDP

END

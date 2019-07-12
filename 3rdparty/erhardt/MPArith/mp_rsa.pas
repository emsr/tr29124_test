unit mp_rsa;

{MP functions for basic RSA based public-key cryptography}

interface

{$i STD.INC}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  MP functions for basic RSA based public-key cryptography

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REMARK        :  In order to keep MPArith (relative) small and modular the
                  RSA sign/verify operations do not call the actual hash
                  functions. The user has to supply a hash algorithm ID
                  and externally calculated hash digests.

 REFERENCES    :  [1] LibTomMath 0.30+ by Tom St Denis
                  [2] MPI by M.J. Fromberger
                   -  RFC 2313 - PKCS #1: RSA Encryption Version 1.5 and
                   -  RFC 3447 - PKCS #1: RSA Encryption Version 2.1 available
                      online from http://tools.ietf.org/html/rfc2313 and
                      http://tools.ietf.org/html/rfc3447
                 [30] J. v. zur Gathen, J. Gerhard, Modern computer algebra, 2nd ed., 2003
                      https://cosec.bit.uni-bonn.de/science/mca/

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.0.01   10.11.05  W.Ehrhardt  Initial version: working pkcs1v15_decrypt a la MPI
 0.0.02   12.11.05  we          EME-PKCS1-v1_5 decoding with $00, $02 ...
 0.0.03   12.11.05  we          Working pkcs1v15_encrypt
 0.0.04   13.11.05  we          nonzero random padding, more error checking
 0.0.05   13.11.05  we          max lengths in more uniform parameter lists
 0.0.06   19.11.05  we          mp_rsa_keygen1, RSA_MINSIZE
 0.0.07   20.11.05  we          mp_rsa_keygen1: changed pt_3mod4 to pt_normal
 0.0.08   11.08.06  we          Avoid FPC warnings: mp_pkcs1v15_encrypt, mp_pkcs1v15_decrypt
 0.0.09   12.11.07  we          Fix memory leak(s) if MPC_HaltOnError is not defined

 1.6.00   22.05.08  we          TPrivateKey, mp_rsa_keygen2, mp_rsa_init/clear/calc_private
 1.6.01   23.05.08  we          mp_rsadp2
 1.6.02   25.05.08  we          mp_i2pchar
 1.6.03   01.06.08  we          mp_rsa_calc_npq, mp_pkcs1v15_decrypt2, mp_rsa_calc_d

 1.7.00   30.06.08  we          correct exception strings in mp_rsa_calc_npq

 1.9.00   29.11.08  we          mp_rsa_calc_nd, fix silly "fillchar(buf^" bug
 1.9.01   02.12.08  we          Uses BTypes: char8, pchar8
 1.9.02   13.12.08  we          mp_rsa_recover_pq
 1.9.03   23.12.08  we          mp_rsa_recover_pq uses s_mp_mca_alg1816
 1.9.04   06.01.09  we          Uses BTypes moved to implementation

 1.10.00  21.01.09  we          changes related to (s)mp_divrem
 1.10.01  10.02.09  we          mp_rsa_recover_pq2

 1.12.00  20.06.09  we          updated RFC URL(s)
 1.12.01  28.06.09  we          improved mp_rsa_calc_npq
 1.12.03  29.07.09  we          Increased trace level in mp_rsa_calc_npq

 1.13.00  10.08.09  we          mp_rsa_wiener

 1.16.00  06.06.10  we          ensure |p-q| is not too small in mp_rsa_calc_npq

 1.17.00  27.12.10  we          Sign/verify functions
 1.17.01  31.12.10  we          mp_pkcs1v15_emsa_encode: SHA224 and RIPEMD-160
 1.17.02  01.01.11  we          Fix wrong byte in RFC4880 SHA224 algorithm identifier

 1.23.00  24.09.12  we          s_mp_mca_alg1816 from mp_numth

 1.33.00  09.06.15  we          moved some mp_clear[x] into success blocks

 1.34.00  08.05.17  we          update some broken links

**************************************************************************)

(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - LibTomMath 0.30+ by Tom St Denis
   - MPI 1.8.6 by Michael J. Fromberger
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*---------------------------------------------------------------------------
Acknowledgement
  RFCs 2313/3447 are based on contributions of RSA Laboratories.
  RSA Security Inc. requests that the PKCS references are identified as
  "RSA Security Inc. PKCS #1 v1.5" and "RSA Security Inc. PKCS #1 v2.1".
----------------------------------------------------------------------------*)

(*-------------------------------------------------------------------------
 (C) Copyright 2005-2017 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)

const
  RSA_MINSIZE = 11; {minimum octet size of modulus}


type
  TPrivateKey = record {RSA private key record with CRT coefficients}
                  p,q  : mp_int; {primes with n=p*q  }
                  dp   : mp_int; {dp = e^-1 mod (p-1)}
                  dq   : mp_int; {dq = e^-1 mod (q-1)}
                  qinv : mp_int; {qinv = q^-1 mod p  }
                end;

type
  TSSAHash = (SA_MD2, SA_MD5, SA_RMD160, SA_SHA1, SA_SHA224, SA_SHA256, SA_SHA384, SA_SHA512);
             {Valid hash algorithms for RSA signature schemes.  Only the hash}
             {digests and algorithm identifiers are used, not the actual hash}
             {functions, i.e. the digests have to be calculated externally.  }


procedure mp_i2osp(const x: mp_int; outp: pointer; len: word);
  {-Convert a nonnegative mp_int to an octet string of a specified length}

procedure mp_i2pchar(const x: mp_int; outp: pointer; len: word);
  {-Convert a nonnegative mp_int to a pchar with a specified max length}

procedure mp_os2ip(inp: pointer; ilen: word; var x: mp_int);
  {-Convert an octet string of length ilen to a nonnegative mp_int}

function  mp_pkcs1v15_decode(em, msg: pointer; emlen: word; var mlen: word): boolean;
  {-EME-PKCS1-v1_5 decoding; true if decoding is successful, false otherwise}

procedure mp_pkcs1v15_decrypt(const d, n: mp_int; ctp,ptp: pointer; clen,pmax: word; var plen: word);
  {-Decrypt a message using RSA and EME-PKCS1-v1_5 padding}

procedure mp_pkcs1v15_decrypt2(const prk: TPrivateKey; const n: mp_int; ctp,ptp: pointer; clen,pmax: word; var plen: word);
  {-Decrypt a message using RSA/CRT and EME-PKCS1-v1_5 padding}

function  mp_pkcs1v15_emsa_encode(AHash: TSSAHash; hdp, smp: pointer; hlen, slen: word): boolean;
  {-EMSA-PKCS1-v1_5 encoding; true if encoding is successful, false otherwise}

function  mp_pkcs1v15_encode(msg, em: pointer; mlen, emlen: word; rnd: byte): boolean;
  {-EME-PKCS1-v1_5 encoding; true if encoding is successful, false otherwise}

procedure mp_pkcs1v15_encrypt(const e, n: mp_int; rnd: byte; ptp,ctp: pointer; plen,cmax: word; var clen: word);
  {-Encrypt a message using RSA and EME-PKCS1-v1_5 padding}

function  mp_pkcs1v15_maxlen(const n: mp_int): word;
  {-Maximum message length for RSA modulus n using PKCS1-v1_5 encoding}

procedure mp_pkcs1v15_sign(const d, n: mp_int; AHash: TSSAHash; hdp,smp: pointer; hlen,smax: word; var slen: word);
  {-Sign a hash digest using RSA and EMSA-PKCS1-v1_5 encoding}

procedure mp_pkcs1v15_sign2(const prk: TPrivateKey; const n: mp_int; AHash: TSSAHash;
                            hdp,smp: pointer; hlen,smax: word; var slen: word);
  {-Sign a hash digest using RSA/CRT and EMSA-PKCS1-v1_5 encoding}

function  mp_pkcs1v15_verify(const e, n: mp_int; AHash: TSSAHash; hdp,smp: pointer; hlen,slen: word): boolean;
  {-Signature verification operation}

procedure mp_rsadp(const c, d, n: mp_int; var m: mp_int);
  {-Basic RSA decryption operation, m=c^d mod n.}

procedure mp_rsadp2(const c: mp_int; const prk: TPrivateKey; var m: mp_int);
  {-Basic RSA decryption operation for private key CRT record.}

procedure mp_rsaep(const m, e, n: mp_int; var c: mp_int);
  {-Basic RSA encryption operation, c=m^e mod n.}

procedure mp_rsasp(const m, d, n: mp_int; var s: mp_int);
  {-Basic RSA signature primitive, s=m^d mod n.}

procedure mp_rsasp2(const m: mp_int; const prk: TPrivateKey; var s: mp_int);
  {-Basic RSA signature primitive for private key CRT record.}

procedure mp_rsavp(const s, e, n: mp_int; var m: mp_int);
  {-Basic RSA verification operation, m=s^e mod n.}

procedure mp_rsa_calc_d(const e: mp_int; const prk: TPrivateKey; var d: mp_int);
  {-Get RSA decryption exponent d from private key record and encryption exponent e}

procedure mp_rsa_calc_nd(const e,p,q: mp_int; var n,d: mp_int);
  {-Calculate n,d from e,p,q}

procedure mp_rsa_calc_npq(const e: mp_int; osize: word; var n,p,q: mp_int);
  {-Generate RSA primes p,q; q<p, n=p*q, osize: octet size of n;}
  { e: public key encryption exponent, e odd and greater than 2}

procedure mp_rsa_calc_private(const e, p, q: mp_int; var prk: TPrivateKey);
  {-Calculate remaining fields of private RSA/CRT key from e,p,q}

procedure mp_rsa_clear_private(var prk: TPrivateKey);
  {-Clear fields of private RSA/CRT key}

procedure mp_rsa_init_private(var prk: TPrivateKey);
  {-Initialize fields of private RSA/CRT key}

procedure mp_rsa_keygen1(const e: mp_int; osize: word; var d, n: mp_int);
  {-Basic RSA private key pair (n, d) generation; osize: octet size of n;}
  { e: public key encryption exponent, e odd and greater than 2}

procedure mp_rsa_keygen2(const e: mp_int; osize: word; var n: mp_int; var prk: TPrivateKey);
  {-Generate private RSA/CRT key prk and modulus n; osize: octet size of n;}
  { e: public key encryption exponent, e odd and greater than 2}

procedure mp_rsa_recover_pq(const n,e,d: mp_int; var p,q: mp_int; var fail: boolean);
  {-Try to recover p,q from n,e,d. Assumes n=p*q with odd primes p,q and}
  { e*d=1 mod lcm(p-1,q-1). Fail=true if no success or if e*d is even.}

procedure mp_rsa_recover_pq2(const n,e,dp: mp_int; var p,q: mp_int; var fail: boolean);
  {-Try to recover p,q from n,e,dp (dp is CRT exponent of p). Assumes n=p*q}
  { with odd primes p,q. Fail=true if no success or if e*dp is even.}

procedure mp_rsa_wiener(const e,n: mp_int; var p,q,d: mp_int; var fail: boolean);
  {-Wiener's attack on small RSA secret exponents: Recover p,q,d from e,n.}
  { Assumes n=p*q with odd primes p>q, e*d=1 mod lcm(p-1,q-1). Fail=true if}
  { no success: typically if bitsize(d) > bitsize(n)/4. If e is given and}
  { d is calculated, usually d will be 'large' and cannot be recovered. }


procedure s_mp_mca_alg1816(const n,L: mp_int; k:integer; var p: mp_int; var fail: boolean);
  {-Try to find a factor p of a squarefree odd integer n>5, with L a multiple}
  { of lambda(n) [lambda: Carmichael function] and a confidence parameter k.}
  { Fail=true, if no success (e.g. n is prime), n is even, n<=5, or L is odd.}


implementation


uses
  BTypes, mp_base, mp_prng, mp_modul, mp_numth;

type
  TByteArr = packed array[0..65519] of byte;  {Helper types}
  PByteArr = ^TByteArr;


{---------------------------------------------------------------------------}
function mp_pkcs1v15_maxlen(const n: mp_int): word;
  {-Maximum message length for RSA modulus n using PKCS1-v1_5 encoding}
var
  modlen: word;
begin
  {MPC_ArgCheck in mp_unsigned_bin_size or deeper}
  modlen := mp_unsigned_bin_size(n);
  {at least RSA_MINSIZE bytes are required for EME-PKCS1-v1_5 padding}
  if modlen<RSA_MINSIZE then mp_pkcs1v15_maxlen := 0
  else mp_pkcs1v15_maxlen := modlen-RSA_MINSIZE;
end;


{---------------------------------------------------------------------------}
procedure mp_i2osp(const x: mp_int; outp: pointer; len: word);
  {-Convert a nonnegative mp_int to an octet string of a specified length}
var
  lx,l0: word;
begin
  {MPC_ArgCheck in mp_unsigned_bin_size or deeper}
  lx := mp_unsigned_bin_size(x);
  if lx>len then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_i2osp: x too large"');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  l0 := len-lx;
  {zero "leading digits"}
  fillchar(outp^,l0,0);
  inc(Ptr2Inc(outp),l0);
  {$ifopt X-} l0 := {$endif} mp_to_unsigned_bin_n(x,outp^,lx);
end;


{---------------------------------------------------------------------------}
procedure mp_i2pchar(const x: mp_int; outp: pointer; len: word);
  {-Convert a nonnegative mp_int to a pchar with a specified max length}
var
  lx,l0: word;
begin
  {MPC_ArgCheck in mp_unsigned_bin_size or deeper}
  lx := mp_unsigned_bin_size(x);
  if lx+1>len then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_i2pchar: x too large"');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  l0 := mp_to_unsigned_bin_n(x,outp^,lx);
  if l0<len then begin
    inc(Ptr2Inc(outp),l0);
    pchar8(outp)^ := #0;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_os2ip(inp: pointer; ilen: word; var x: mp_int);
  {-Convert an octet string of length ilen to a nonnegative mp_int}
begin
  {MPC_ArgCheck in mp_read_unsigned_bin}
  mp_read_unsigned_bin(x, inp^, ilen);
end;


{---------------------------------------------------------------------------}
function  mp_pkcs1v15_decode(em, msg: pointer; emlen: word; var mlen: word): boolean;
  {-EME-PKCS1-v1_5 decoding; true if decoding is successful, false otherwise}
  { em    - encoded message
    emlen - length of encoded message, in bytes
    msg   - decode message (may be same as em)
    mlen  - length of decoded msg}
var
  ep: PByteArr absolute em;
  ip: word;
  res: boolean;
begin
  {EME-PKCS1-v1_5 decoding: Separate the encoded message EM into an
   octet string PS consisting of nonzero octets and a message M as

      EM = 0x00 || 0x02 || PS || 0x00 || MSG.

   If the first octet of EM does not have hexadecimal value 0x00, if the
   second octet of EM does not have hexadecimal value 0x02, if there is
   no octet with hexadecimal value 0x00 to separate PS from M, or if the
   length of PS is less than 8 octets, output "decryption error".

   Note: Care shall be taken to ensure that an opponent cannot distinguish
   the different error conditions, whether by error message or timing.
   Otherwise an opponent may be able to obtain useful information about
   the decryption of the ciphertext C. WE note: This is partly achieved
   by running through the entire code without error exit.}

  res := true;

  {check min length and starting bytes}
  if emlen<RSA_MINSIZE then res := false;
  if ep^[0]<>$00 then res := false;
  if ep^[1]<>$02 then res := false;

  {look for zero separator}
  ip := 2;
  while ip<emlen do begin
    if ep^[ip]=0 then break;
    inc(ip);
  end;

  {position of zero separator must be >= 10 and < emlen}
  if (ip=emlen) or (ip<10) then res := false;

  {copy msg}
  mlen := 0;
  if emlen > (ip + 1) then begin
    mlen := emlen - (ip + 1);
    move(ep^[ip+1],msg^,mlen);
  end
  else res := false;

  mp_pkcs1v15_decode := res;

end;


{---------------------------------------------------------------------------}
procedure mp_pkcs1v15_decrypt(const d, n: mp_int; ctp,ptp: pointer; clen,pmax: word; var plen: word);
  {-Decrypt a message using RSA and EME-PKCS1-v1_5 padding}
  { n,d   - RSA private key pair (n, d) = (modulus, exponent)
    ctp   - input message (ciphertext)
    ptp   - pointer to buffer holding decrypted plaintext
    clen  - length of input message, in bytes
    pmax  - max. length of plaintext buffer
    plen  - length of plaintext, in bytes}
var
  k: word;
  c: mp_int;
  buf: pointer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    {check d here; n checked in mp_unsigned_bin_size}
    if mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pkcs1v15_decrypt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if clen<>k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_decrypt: clen <> k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  if pmax<k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_decrypt: pmax < k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  {get local memory for encoded message}
  buf := mp_getmem(clen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_decrypt: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {initialize ciphertext mp_int representative}
  mp_init(c);
  if mp_error=MP_OKAY then begin
    {convert ciphertext to mp_int representative}
    mp_os2ip(ctp, clen, c);
    {apply the RSADP decryption primitive to the ciphertext representative c}
    mp_rsadp(c, d, n, c);
    {Convert the message representative to an encoded message}
    mp_i2osp(c, buf, k);
    {EME-PKCS1-v1_5 decoding buf -> ptp}
    if not mp_pkcs1v15_decode(buf, ptp, k, plen) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXUndef.Create('mp_pkcs1v15_decrypt: decode error');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$endif}
      set_mp_error(MP_UNDEF);
    end;
    mp_clear(c);
  end;
  fillchar(buf^,clen,0);
  mp_freemem(buf, clen);
end;


{---------------------------------------------------------------------------}
procedure mp_pkcs1v15_decrypt2(const prk: TPrivateKey; const n: mp_int; ctp,ptp: pointer; clen,pmax: word; var plen: word);
  {-Decrypt a message using RSA/CRT and EME-PKCS1-v1_5 padding}
  { prk,n - private key CRT record and modulus
    ctp   - input message (ciphertext)
    ptp   - pointer to buffer holding decrypted plaintext
    clen  - length of input message, in bytes
    pmax  - max. length of plaintext buffer
    plen  - length of plaintext, in bytes}
var
  k: word;
  c: mp_int;
  buf: pointer;
begin
  if mp_error<>MP_OKAY then exit;
  {Init check n in mp_unsigned_bin_size or below}
  {get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if clen<>k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_decrypt2: clen <> k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  if pmax<k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_decrypt2: pmax < k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  {get local memory for encoded message}
  buf := mp_getmem(clen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_decrypt2: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {initialize ciphertext mp_int representative}
  mp_init(c);
  if mp_error=MP_OKAY then begin
    {convert ciphertext to mp_int representative}
    mp_os2ip(ctp, clen, c);
    {apply the RSADP decryption primitive to the ciphertext representative c}
    mp_rsadp2(c, prk, c);
    {Convert the message representative to an encoded message}
    mp_i2osp(c, buf, k);
    {EME-PKCS1-v1_5 decoding buf -> ptp}
    if not mp_pkcs1v15_decode(buf, ptp, k, plen) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXUndef.Create('mp_pkcs1v15_decrypt2: decode error');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$endif}
      set_mp_error(MP_UNDEF);
    end;
    mp_clear(c);
  end;
  fillchar(buf^,clen,0);
  mp_freemem(buf, clen);
end;


{---------------------------------------------------------------------------}
function mp_pkcs1v15_encode(msg, em: pointer; mlen, emlen: word; rnd: byte): boolean;
  {-EME-PKCS1-v1_5 encoding; true if encoding is successful, false otherwise}
  { msg   - plaintext message
    mlen  - length of plaintext message
    em    - encoded message
    emlen - length of encoded message, in bytes, will be padded
    rnd   - pad byte, if pad=0 use bytes generated with mp_random_byte}
var
  i,ip: word;
  ep: PByteArr absolute em;
  rb: byte;
begin
  if mlen>emlen-RSA_MINSIZE then mp_pkcs1v15_encode := false
  else begin
    {EME-PKCS1-v1_5 encoding:
     Generate an octet string PS of length emLen - mLen - 3 consisting
     of pseudo-randomly generated nonzero octets.  The length of PS
     will be at least eight octets.

     Concatenate PS, the message MSG, and other padding to form an
     encoded message EM of length k octets as

           EM = 0x00 || 0x02 || PS || 0x00 || MSG.
     WE-Note: First 0 byte makes mp_int representative less than modulus}
    {insert starting bytes}
    ep^[0] := $00;
    ep^[1] := $02;
    {calculate padding offset. Note: mlen<=emlen-RSA_MINSIZE -> ip>=RSA_MINSIZE}
    ip := emlen-mlen;
    {insert padding bytes}
    for i:=2 to ip-2 do begin
      if rnd<>0 then ep^[i] := rnd
      else begin
        {generate random nonzero padding bytes}
        repeat
          rb := mp_random_byte;
        until rb<>0;
        ep^[i] := rb;
      end;
    end;
    {insert zero separator}
    ep^[ip-1] := 0;
    move(msg^,ep^[ip],mlen);
    mp_pkcs1v15_encode := true;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_pkcs1v15_encrypt(const e, n: mp_int; rnd: byte; ptp,ctp: pointer; plen,cmax: word; var clen: word);
  {-Encrypt a message using RSA and EME-PKCS1-v1_5 padding}
  { n,e   - recipient's RSA public key pair (n, e) = (modulus, exponent)
    rnd   - padding byte, if 0 mp_random_byte is used
    ptp   - input message (plaintext)
    ctp   - output message (ciphertext)
    plen  - length of buffer holding plaintext, in bytes
    cmax  - max length of ciphertext buffer
    clen  - length of output message, in bytes}
var
  k: word;
  c: mp_int;
  buf: pointer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    {check e here; n checked in mp_unsigned_bin_size}
    if mp_not_init(e) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pkcs1v15_encrypt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if cmax<k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_encrypt: cmax < k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  clen := k;

  {get local memory for encoded message}
  buf := mp_getmem(clen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_encrypt: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {initialize ciphertext mp_int representative}
  mp_init(c);
  if mp_error=MP_OKAY then begin
    {-EME-PKCS1-v1_5 encoding ptp -> buf}
    if not mp_pkcs1v15_encode(ptp, buf, plen, k, rnd) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXUndef.Create('mp_pkcs1v15_encrypt: encode error');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$endif}
      set_mp_error(MP_UNDEF);
    end
    else begin
      {Convert the encoded message EM to mp_int representative c}
      mp_os2ip(buf, k, c);
      {Apply the RSAEP encryption primitive to the message representative c}
      mp_rsaep(c, e, n, c);
      {Convert the ciphertext representative c to octet string}
      mp_i2osp(c, ctp, clen);
    end;
    mp_clear(c);
  end;
  fillchar(buf^,clen,0);
  mp_freemem(buf, clen);
end;


{---------------------------------------------------------------------------}
procedure mp_rsadp(const c, d, n: mp_int; var m: mp_int);
  {-Basic RSA decryption operation, m=c^d mod n.}
  { c   : ciphertext representative, an integer between 0 and n - 1.
    n,d : RSA private key pair (n, d) = (modulus, exponent).
    m   : message representative, an integer between 0 and n - 1.}
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    {check c,n here; d,m are checked in mp_exptmod}
    if mp_not_init(c) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsadp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Insure that ciphertext representative is in range of modulus}
  if (c.sign=MP_NEG) or (mp_cmp(c, n)<>MP_LT) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rsadp: ciphertext representative out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_exptmod(c, d, n, m);
end;


{---------------------------------------------------------------------------}
procedure mp_rsadp2(const c: mp_int; const prk: TPrivateKey; var m: mp_int);
  {-Basic RSA decryption operation for private key CRT record.}
  { c   : ciphertext representative, an integer between 0 and n - 1. Note
          that this condition will be checked only approximately.
    prk : RSA private key CRT record.
    m   : message representative, an integer between 0 and n - 1.}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    {check c,n here; d,m are checked in mp_exptmod}
    if mp_not_init(c) or mp_not_init(m) or mp_not_init(prk.p) or mp_not_init(prk.q)
       or mp_not_init(prk.dp) or mp_not_init(prk.dq) or mp_not_init(prk.qinv)
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsadp2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Insure that ciphertext representative is in range of modulus}
  if (c.sign=MP_NEG) or (mp_bitsize(c) > (mp_bitsize(prk.p)+mp_bitsize(prk.q))) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rsadp2: ciphertext representative out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init(t);
  if mp_error=MP_OKAY then with prk do begin
    {no need to reduce c mod p/q, will be done implicitly in mp_exptmod}
    {m1 = c^dp mod p}
    mp_exptmod(c,dp,p,t);
    {m2 = c^dq mod q}
    mp_exptmod(c,dq,q,m);
    {h  = (m1-m2)*qinv mod p}
    mp_submod(t,m,p,t);
    mp_mulmod(t,qinv,p,t);
    {m = m2+q*h}
    mp_mul(q,t,t);
    mp_add(t,m,m);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsaep(const m, e, n: mp_int; var c: mp_int);
  {-Basic RSA encryption operation, c=m^e mod n.}
  { m   : message representative, an integer between 0 and n - 1.
    n,e : RSA public key pair (n, e) = (modulus, exponent).
    c   : ciphertext representative, an integer between 0 and n - 1.}
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    {check m,n here; d,m are checked in mp_exptmod}
    if mp_not_init(m) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsaep');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Insure that message representative is in range of modulus}
  if (m.sign=MP_NEG) or (mp_cmp(m, n)<>MP_LT) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rsaep: message representative out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_exptmod(m, e, n, c);
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_calc_d(const e: mp_int; const prk: TPrivateKey; var d: mp_int);
  {-Get RSA decryption exponent d from private key record and encryption exponent e}
var
  p,q: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(e) or mp_not_init(d) or mp_not_init(prk.p) or mp_not_init(prk.q)
       or mp_not_init(prk.dp) or mp_not_init(prk.dq) or mp_not_init(prk.qinv)
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_calc_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init2(p,q);
  if mp_error=MP_OKAY then begin
    {now calculate d := e^-1 mod lcm(p-1, q-1)}
    mp_sub_d(prk.p,1,p);
    mp_sub_d(prk.q,1,q);
    mp_lcm(p,q,d);
    mp_invmod(e,d,d);
    mp_clear2(p,q);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_calc_nd(const e,p,q: mp_int; var n,d: mp_int);
  {-Calculate n,d from e,p,q}
var
  p1,q1: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(e) or mp_not_init(d) or mp_not_init(n) or mp_not_init(p) or mp_not_init(q) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_calc_nd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init2(p1,q1);
  if mp_error=MP_OKAY then begin
    mp_sub_d(p,1,p1);
    mp_sub_d(q,1,q1);
    {n := p*q}
    mp_mul(p,q,n);
    {d := e^-1 mod lcm(p-1, q-1)}
    mp_lcm(p1,q1,d);
    mp_invmod(e,d,d);
    mp_clear2(p1,q1);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_calc_npq(const e: mp_int; osize: word; var n,p,q: mp_int);
  {-Generate RSA primes p,q; q<p, n=p*q, osize: octet size of n;}
  { e: public key encryption exponent, e odd and greater than 2}
var
  t: mp_int;
  s3,s4,s8: word;
  f: mp_digit;
const
  fmaxs = mp_digit(MP_DIGIT_MAX and $3FF);
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(e) or mp_not_init(n) or mp_not_init(p) or mp_not_init(q) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_calc_npq');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {check if e is odd and greater than 2}
  if mp_iseven(e) or (mp_cmp_d(e,3)=MP_LT) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_rsa_calc_npq: e is even or <3');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {check requested octet size of n}
  if osize<RSA_MINSIZE then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_rsa_calc_npq: osize < RSA_MINSIZE');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  mp_init(t);
  if mp_error=MP_OKAY then begin
    {bit size of p and q = 1/2 bit size of n = 4*osize}
    s3 := 3*osize;
    s4 := 4*osize;
    s8 := s4+s4;
    {generate p with gcd(p-1,e)=1}
    repeat
      if mp_error<>MP_OKAY then begin
        mp_clear(t);
        exit;
      end;
      mp_rand_bits(p,s4);
      p.pdigits^[0] := p.pdigits^[0] or 1;
      {first test if p has small factor}
      mp_small_factor(p,3,fmaxs,f);
      if f<>0 then continue;
      {now check gcd(p-1,e)=1}
      mp_sub_d(p,1,t);
      if mp_gcd1(e,t,t) then begin
        {done if p is a BPSW probable prime}
        if mp_is_spsp_d(p,2) and mp_is_slpsp(p) then break;
      end;
    until false;
    {generate q with gcd(q-1,e)=1}
    repeat
      if mp_error<>MP_OKAY then begin
        mp_clear(t);
        exit;
      end;
      mp_rand_bits(q,s4);
      q.pdigits^[0] := q.pdigits^[0] or 1;
      {ensure |p-q| is not too small}
      mp_sub(p,q,t);
      if mp_bitsize(t) < s3 then begin
        continue;
      end;
      {test if q has a small factor}
      mp_small_factor(q,3,fmaxs,f);
      if f<>0 then continue;
      {check if n=p*q has required bit size}
      mp_mul(p,q,n);
      if mp_bitsize(n)<>s8 then continue;
      {now check gcd(q-1,e)=1}
      mp_sub_d(q,1,t);
      if mp_gcd1(e,t,t) then begin
        {done if q is a BPSW probable prime}
        if mp_is_spsp_d(q,2) and mp_is_slpsp(q) then break;
      end;
    until false;
    if mp_is_le(p,q) then mp_exch(p,q);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_calc_private(const e, p, q: mp_int; var prk: TPrivateKey);
  {-Calculate remaining fields of private RSA/CRT key from e,p,q}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(e) or mp_not_init(p) or mp_not_init(q) or mp_not_init(prk.p) or mp_not_init(prk.q)
       or mp_not_init(prk.dp) or mp_not_init(prk.dq) or mp_not_init(prk.qinv)
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_calc_private');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_copy(p, prk.p);
  mp_copy(q, prk.q);
  with prk do begin
    {dp := e^-1 mod (p-1)}
    mp_sub_d(p,1,qinv);
    mp_invmod(e,qinv,dp);
    {dq := e^-1 mod (q-1)}
    mp_sub_d(q,1,qinv);
    mp_invmod(e,qinv,dq);
    {qinv := q^-1 mod p}
    mp_invmod(q,p,qinv);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_clear_private(var prk: TPrivateKey);
  {-Clear fields of private RSA/CRT key}
begin
  with prk do mp_clear5(p,q,dp,dq,qinv);
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_init_private(var prk: TPrivateKey);
  {-Initialize fields of private RSA/CRT key}
begin
  with prk do mp_init5(p,q,dp,dq,qinv);
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_keygen1(const e: mp_int; osize: word; var d, n: mp_int);
  {-Basic RSA private key pair (n, d) generation; osize: octet size of n}
  { e: public key encryption exponent, e odd and greater than 2}
var
  p,q: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  mp_init2(p,q);
  if mp_error=MP_OKAY then begin
    {first get RSA primes p,q}
    mp_rsa_calc_npq(e,osize,n,p,q);
    {now calculate d := e^-1 mod lcm(p-1, q-1)}
    mp_dec(p);
    mp_dec(q);
    mp_lcm(p,q,d);
    mp_invmod(e,d,d);
    mp_clear2(p,q);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_keygen2(const e: mp_int; osize: word; var n: mp_int; var prk: TPrivateKey);
  {-Generate private RSA/CRT key prk and modulus n; osize: octet size of n;}
  { e: public key encryption exponent, e odd and greater than 2}
begin
  if mp_error<>MP_OKAY then exit;
  with prk do begin
    {$ifdef MPC_ArgCheck}
      {other args ar checked in mp_rsa_calc_npq}
      if mp_not_init(dp) or mp_not_init(dq) or mp_not_init(qinv) then begin
        {$ifdef MPC_HaltOnArgCheck}
          {$ifdef MPC_UseExceptions}
            raise MPXNotInit.Create('mp_rsa_keygen2');
          {$else}
            RunError(MP_RTE_NOTINIT);
          {$endif}
        {$else}
          set_mp_error(MP_NOTINIT);
          exit;
        {$endif}
      end;
    {$endif}

    {first get RSA primes p,q}
    mp_rsa_calc_npq(e,osize,n,p,q);
    {dp := e^-1 mod (p-1)}
    mp_sub_d(p,1,qinv);
    mp_invmod(e,qinv,dp);
    {dq := e^-1 mod (q-1)}
    mp_sub_d(q,1,qinv);
    mp_invmod(e,qinv,dq);
    {qinv := q^-1 mod p}
    mp_invmod(q,p,qinv);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_recover_pq(const n,e,d: mp_int; var p,q: mp_int; var fail: boolean);
  {-Try to recover p,q from n,e,d. Assumes n=p*q with odd primes p,q and}
  { e*d=1 mod lcm(p-1,q-1). Fail=true if no success or if e*d is even.}
var
  L: mp_int;
begin
  fail := true;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(e) or mp_not_init(d) or mp_not_init(p) or mp_not_init(q) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_recover_pq');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init(L);
  if mp_error=MP_OKAY then begin
    {L=e*d-1 is a multiple of lambda(n)=lcm(p-1,q-1)}
    mp_mul(d,e,L);
    mp_dec(L);
    s_mp_mca_alg1816(n,L,100,p,fail);
    if not fail then begin
      mp_div(n,p,q);
      {normalize p>q}
      if mp_cmp_mag(p,q)<>MP_GT then mp_exch(p,q);
      fail := mp_error<>MP_OKAY;
    end;
    mp_clear(L);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_mca_alg1816(const n,L: mp_int; k:integer; var p: mp_int; var fail: boolean);
  {-Try to find a factor p of a squarefree odd integer n>5, with L a multiple}
  { of lambda(n) [lambda: Carmichael function] and a confidence parameter k.}
  { Fail=true, if no success (e.g. n is prime), n is even, n<=5, or L is odd.}
label
  leave;
var
  x,y: mp_int;
  h: integer;
  i,t: longint;
begin
  {Uses Algorithm 18.16 (Special integer factorization) from the}
  {solution to exercise 18.12 (ii) of MCA [30], available online}
  {as https://cosec.bit.uni-bonn.de/fileadmin/user_upload/science/mca/solutions.pdf}

  fail := true;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(L) or mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_mca_alg1816');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Easy outs}
  if mp_isodd(L) or mp_iseven(n) or (mp_cmp_d(n,5)<0) then exit;

  mp_init2(x,y);
  if mp_error<>MP_OKAY then exit;

  {The variables x,y are used for varying entities from Alg. 18.16; x will}
  {be a non-trivial factor found by gcd(n,.) if fail=false}

  for h:=1 to k do begin
    if mp_error<>MP_OKAY then goto leave;
    {choose 1 < a < n-1 uniformly at random}
    mp_sub_d(n,3,x);   {OK, because n>6}
    mp_rand(y,n.used);
    mp_mod(y,x,y);     {0 <= y < n-3}
    mp_add_d(y,2,y);   {2 <= y < n-1}

    {check if a factor is found by accident}
    if not mp_gcd1(y,n,x) then begin
      {x<>1. Since y < n, x will be a proper factor}
      fail := false;
      goto leave;
    end;

    {write L = 2^t*m with m odd and t>0}
    mp_makeodd(L,x,t);

    {calculate b0 = a^m mod n}
    mp_exptmod(y, x, n, y);

    {if b0=1, try another random a}
    if mp_is1(y) then continue;

    {here y is the last b_i with b_i <> 1}
    for i:=1 to t do begin
      {calculate x=b_i}
      mp_sqrmod(y,n,x);
      if mp_is1(x) then begin
        {x=b_i=1, i.e. found maximal j=i-1 with b_j<>1, use it for trial gcd}
        mp_inc(y);
        {exclude trivial case b_j+1=n}
        if mp_is_ne(n,y) then begin
          if not mp_gcd1(y,n,x) then begin
            fail := false;
            goto leave;
          end;
        end;
        {skip useless squarings because all other b_i will be = 1}
        break; {for i loop}
      end;
      {update b_i with b_i <> 1}
      mp_exch(x,y);
      if mp_error<>MP_OKAY then goto leave;
    end;
  end;

leave:

  if (not fail) and (mp_error=MP_OKAY) then begin
    {x is a proper factor}
    mp_exch(x,p);
  end;

  mp_clear2(x,y);
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_recover_pq2(const n,e,dp: mp_int; var p,q: mp_int; var fail: boolean);
  {-Try to recover p,q from n,e,dp (dp is CRT exponent of p). Assumes n=p*q}
  { with odd primes p,q. Fail=true if no success or if e*dp is even.}
label
  leave;
var
  r,v,x: mp_int;
  w: mp_int absolute v;
  k: integer;
  i,s: longint;
begin

  {This procedure implements the probabilistic polynomial time Algorithm 2}
  {from S. Maitra, S. Sarkar: Polynomial-Time Equivalence of Computing the}
  {CRT-RSA Secret Key(s) and Factoring, http://eprint.iacr.org/2009/062   }

  fail := true;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(e) or mp_not_init(dp) or mp_not_init(p) or mp_not_init(q) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_recover_pq2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_iseven(n) or (mp_cmp_d(n,5)<0) then exit;

  mp_init3(r,v,x);
  if mp_error<>MP_OKAY then exit;

  {Step 1: calculate r,s with e*dp-1=2^s*r}
  mp_mul(dp,e,x);
  mp_dec(x);
  mp_makeodd(x,r,s);

  {if s=0 then p is not odd, because e*dp-1 = m*(p-1) must be even}
  if s=0 then goto leave;

  {Algorithm 2 is repeated max k=100 times}
  for k:=1 to 100 do begin
    {generate w in the range 2 <= w < n}
    repeat
      mp_rand(w, n.used);
      mp_mod(w,n,v);
      if mp_error<>MP_OKAY then goto leave;
    until mp_cmp_d(w,1)=MP_GT;
    {Steps 2,3}
    if not mp_gcd1(w,n,x) then begin
      fail := false;
      goto leave;
    end;
    {Step 4}
    for i:=0 to s do begin
      {4a: v := w^(2^i*r)}
      if i=0 then mp_exptmod(w,r,n,v)
      else mp_sqrmod(v,n,v);
      {4b: if v=0 or v=1 then failure, i.e. try next random w}
      mp_sub_d(v,1,x);
      if mp_is0(v) or mp_is0(x) then break;
      {4c, 4d: calculate gcd(v-1,n), success if 1 < gcd < n}
      if not mp_gcd1(x,n,x) then begin
        fail := false;
        goto leave;
      end;
      if mp_error<>MP_OKAY then goto leave;
    end;
  end;

leave:

  if (not fail) and (mp_error=MP_OKAY) then begin
    {x is non trivial factor, set q=n/x, p=x}
    mp_div(n,x,q);
    mp_exch(x,p);
    if mp_is1(p) or mp_is1(q) then fail := true;
  end;

  mp_clear3(r,v,x);
end;


{---------------------------------------------------------------------------}
procedure mp_rsa_wiener(const e,n: mp_int; var p,q,d: mp_int; var fail: boolean);
  {-Wiener's attack on small RSA secret exponents: Recover p,q,d from e,n.}
  { Assumes n=p*q with odd primes p>q, e*d=1 mod lcm(p-1,q-1). Fail=true if}
  { no success: typically if bitsize(d) > bitsize(n)/4. If e is given and}
  { d is calculated, usually d will be 'large' and cannot be recovered. }
var
  r1,r2,a1,a2,b1,b2,t,k,dg: mp_int;
  i: integer;
begin
  { M. J. Wiener, Cryptanalysis of short RSA secret exponents, 1990, IEEE}
  { Transactions on Information Theory, vol. 36, p.553-558, available via}
  { http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.92.5261}

  {e*d = 1 mod lcm(p-1,q-1) implies e*d = 1 + K*lcm(p-1,q-1). Define}
  {G,g,k as G = gcd(p-1,q-1), g=G/gcd(K,G), k=K/gcd(K,G). It follows}
  {e*d = 1 + (K/G)*(p-1)*(q-1) = 1 + (k/g)*(p-1)*(q-1). Divide by dpq:}

  {e/n = e/(p*q) = k/(dgpq)*(p-1)*(q-1) + 1/(dpg) = k/(dg)*(1+x)}
  {where  x = (p + q - 1 - g/k) / n is typically 'small'}

  fail := true;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(e) or mp_not_init(d) or mp_not_init(p) or mp_not_init(q) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rsa_wiener');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init9(r1,r2,a1,a2,b1,b2,t,k,dg);

  {Generate continued fraction convergents a(i)/b(i) for e/n}
  {Initialize a1,a2,b1,b2 for the iteration, ie for i=-1,-2}
  mp_zero(a1);
  mp_set1(b1);
  mp_zero(b2);
  mp_set1(a2);

  i:=0;

  {run Euclid's algorithm for r1=|e|, r2=|n|}
  mp_abs(e,r1);
  mp_abs(n,r2);
  while r2.used > 0 do begin
    {q_i = r2 div r1;  k is used to hold q_i}
    {(r2,r1) = (r1, r2 - q_i*r1) = (r1, r2 mod r1)}
    mp_exch(r2,r1);
    mp_divrem(r2,r1,@k,@r2);

    {a(i) = q_i*a(i-1) + a(i-2)}
    mp_exch(a2,a1);
    mp_mul(a1,k,t);
    mp_add(t,a2,a2);

    {b(i) = q_i*b(i-1) + b(i-2)}
    mp_exch(b2,b1);
    mp_mul(b1,k,t);
    mp_add(t,b2,b2);

    {Next guesses for k and dg. If i is even, use q_i+1 in cf convergent}
    {because the guessed quotient should always be >= the true quotient}
    if odd(i) then begin
      mp_copy(a2,k);
      mp_copy(b2,dg);
    end
    else begin
      {a2' = a2 + a1}
      {b2' = b2 + b1}
      mp_add(a2,a1,k);
      mp_add(b2,b1,dg);
    end;

    {if k=0 skip to next iteration}
    if k.used> 0 then begin
      {guess edg = e*dg}
      mp_mul(dg,e,t);
      {guess phi = (p-1)*(q-1) = edg div k and g = edg mod k}
      mp_divrem(t,k,@p,@k);
      {if guessed phi=0 skip to next iteration}
      if p.used>0 then begin
        {from here to end of loop k == g}
        {guess p+q = pq - (p-1)*(q-1) + 1 = n - phi + 1}
        mp_sub(n,p,p);
        mp_inc(p);
        {if p+q is odd skip to next iteration}
        if mp_iseven(p) then begin
          {guess t = ((p-q)/2)^2 = ((p+q)/2)^2 - pq =  ((p+q)/2)^2 - n}
          mp_shr1(p);
          mp_sqr(p,t);
          mp_sub(t,n,t);
          {if t is no perfect square skip to next iteration}
          if mp_is_square2(t,@t) then begin
            {q = (p+q)/2 - (p-q)/2}
            {p = (p+q)/2 + (p-q)/2}
            mp_sub(p,t,q);
            mp_add(p,t,p);
            {d = dg div g}
            mp_divrem(dg,k,@d,@t);
            fail := false;
            {$ifdef MPC_USE_Assert}
              {only if assert supported by compiler or debug}
              assert(mp_is0(t), MPAF+'dg mod g = 0');
            {$endif}
            break;
          end;
        end;
      end;
    end;
    inc(i);
  end;
  mp_clear9(r1,r2,a1,a2,b1,b2,t,k,dg);
end;

{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
procedure mp_rsasp(const m, d, n: mp_int; var s: mp_int);
  {-Basic RSA signature primitive, s=m^d mod n.}
  { m   : message representative, an integer between 0 and n - 1.
    n,d : RSA private key pair (n, d) = (modulus, exponent).
    s   : signature representative, an integer between 0 and n - 1.}
begin
  mp_rsadp(m,d,n,s);
end;


{---------------------------------------------------------------------------}
procedure mp_rsasp2(const m: mp_int; const prk: TPrivateKey; var s: mp_int);
  {-Basic RSA signature primitive for private key CRT record.}
  { m   : message representative, an integer between 0 and n - 1.
    prk : RSA private key CRT record.
    s   : signature representative, an integer between 0 and n - 1. Note
          that this condition will be checked only approximately.}
begin
  mp_rsadp2(m,prk,s)
end;


{---------------------------------------------------------------------------}
procedure mp_rsavp(const s, e, n: mp_int; var m: mp_int);
  {-Basic RSA verification operation, m=s^e mod n.}
  {  s : signature representative, an integer between 0 and n - 1.
     m : message representative, an integer between 0 and n - 1.
    n,e: RSA public key pair (n, e) = (modulus, exponent).}
begin
  mp_rsaep(s,e,n,m);
end;


{---------------------------------------------------------------------------}
function mp_pkcs1v15_emsa_encode(AHash: TSSAHash; hdp, smp: pointer; hlen, slen: word): boolean;
  {-EMSA-PKCS1-v1_5 encoding; true if encoding is successful, false otherwise}
  { AHash - hash algorithm for signature
    hdp   - hash digest
    hlen  - length of hash digest
    smp   - encoded signature message
    slen  - length of encoded signature, in bytes, will be padded}
var
  i,ip,tlen: word;
  ep: PByteArr absolute smp;
  aip: pointer;
(* From RFC 3447, Section 9.2 EMSA-PKCS1-v1_5, Note 1
   MD2:     (0x)30 20 30 0c 06 08 2a 86 48 86 f7 0d 02 02 05 00 04 10 || H.
   MD5:     (0x)30 20 30 0c 06 08 2a 86 48 86 f7 0d 02 05 05 00 04 10 || H.
   SHA-1:   (0x)30 21 30 09 06 05 2b 0e 03 02 1a 05 00 04 14 || H.
   SHA-256: (0x)30 31 30 0d 06 09 60 86 48 01 65 03 04 02 01 05 00 04 20 || H.
   SHA-384: (0x)30 41 30 0d 06 09 60 86 48 01 65 03 04 02 02 05 00 04 30 || H.
   SHA-512: (0x)30 51 30 0d 06 09 60 86 48 01 65 03 04 02 03 05 00 04 40 || H.

 * From RFC 4880, Sec 5.2.2.  Version 3 Signature Packet Format
   (Second byte of SHA224 from Errata ID: 2270, wrong original reads $31)
   MD5:        0x30, 0x20, 0x30, 0x0C, 0x06, 0x08, 0x2A, 0x86,
               0x48, 0x86, 0xF7, 0x0D, 0x02, 0x05, 0x05, 0x00,
               0x04, 0x10
   RIPEMD-160: 0x30, 0x21, 0x30, 0x09, 0x06, 0x05, 0x2B, 0x24,
               0x03, 0x02, 0x01, 0x05, 0x00, 0x04, 0x14
   SHA-1:      0x30, 0x21, 0x30, 0x09, 0x06, 0x05, 0x2b, 0x0E,
               0x03, 0x02, 0x1A, 0x05, 0x00, 0x04, 0x14
   SHA224:     0x30, 0x2d, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86,
               0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x04, 0x05,
               0x00, 0x04, 0x1C
   SHA256:     0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86,
               0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x01, 0x05,
               0x00, 0x04, 0x20
   SHA384:     0x30, 0x41, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86,
               0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x02, 0x05,
               0x00, 0x04, 0x30
   SHA512:     0x30, 0x51, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86,
               0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x03, 0x05,
               0x00, 0x04, 0x40
*)
const
  ai_md2:    array[0..17] of byte = ($30,$20,$30,$0c,$06,$08,$2a,$86,$48,$86,$f7,$0d,$02,$02,$05,$00,$04,$10);
  ai_md5:    array[0..17] of byte = ($30,$20,$30,$0c,$06,$08,$2a,$86,$48,$86,$f7,$0d,$02,$05,$05,$00,$04,$10);
  ai_rmd160: array[0..14] of byte = ($30,$21,$30,$09,$06,$05,$2B,$24,$03,$02,$01,$05,$00,$04,$14);
  ai_sha1:   array[0..14] of byte = ($30,$21,$30,$09,$06,$05,$2b,$0e,$03,$02,$1a,$05,$00,$04,$14);
  ai_sha224: array[0..18] of byte = ($30,$2d,$30,$0d,$06,$09,$60,$86,$48,$01,$65,$03,$04,$02,$04,$05,$00,$04,$1C);
  ai_sha256: array[0..18] of byte = ($30,$31,$30,$0d,$06,$09,$60,$86,$48,$01,$65,$03,$04,$02,$01,$05,$00,$04,$20);
  ai_sha384: array[0..18] of byte = ($30,$41,$30,$0d,$06,$09,$60,$86,$48,$01,$65,$03,$04,$02,$02,$05,$00,$04,$30);
  ai_sha512: array[0..18] of byte = ($30,$51,$30,$0d,$06,$09,$60,$86,$48,$01,$65,$03,$04,$02,$03,$05,$00,$04,$40);
begin

  (* EMSA-PKCS1-v1_5 encoding:
     Generate an octet string PS of length slen-(tlen+hlen) - 3 consisting
     of $FF bytes. The length of PS will be at least eight octets.

     Concatenate PS, the DigestInfo (algorithm identifier || hash digest) and
     other padding to form an signature message SM of length slen octets:
     SM = 0x00 || 0x01 || PS || 0x00 || DigestInfo.

     WE-Note: First 0 byte makes mp_int representative less than modulus
  *)

  mp_pkcs1v15_emsa_encode := false;

  case AHash of
       SA_MD2: begin tlen := sizeof(ai_md2);    aip := @ai_md2;    end;
       SA_MD5: begin tlen := sizeof(ai_md5);    aip := @ai_md5;    end;
    SA_RMD160: begin tlen := sizeof(ai_rmd160); aip := @ai_rmd160; end;
      SA_SHA1: begin tlen := sizeof(ai_sha1);   aip := @ai_sha1;   end;
    SA_SHA224: begin tlen := sizeof(ai_sha224); aip := @ai_sha224; end;
    SA_SHA256: begin tlen := sizeof(ai_sha256); aip := @ai_sha256; end;
    SA_SHA384: begin tlen := sizeof(ai_sha384); aip := @ai_sha384; end;
    SA_SHA512: begin tlen := sizeof(ai_sha512); aip := @ai_sha512; end;
         else  exit; {invalid AHash}
  end;

  {check hash length}
  if hlen<>PByteArr(aip)^[pred(tlen)] then exit;
  {check if slen is large enough}
  if tlen + hlen + RSA_MINSIZE > slen then exit;

  {Here AHash and lengths are valid; calculate padding offset.}
  {Note: hlen + tlen <= slen - RSA_MINSIZE -> ip >= RSA_MINSIZE}
  ip := slen - tlen - hlen;

  {insert starting bytes}
  ep^[0] := $00;
  ep^[1] := $01;
  {insert $FF padding bytes}
  for i:=2 to ip-2 do ep^[i] := $FF;
  {insert zero separator}
  ep^[ip-1] := 0;
  {insert algorithm identifier}
  move(aip^,ep^[ip],tlen);
  {insert hash digest}
  move(hdp^,ep^[ip+tlen],hlen);

  mp_pkcs1v15_emsa_encode := true;
end;


{---------------------------------------------------------------------------}
function mp_pkcs1v15_verify(const e, n: mp_int; AHash: TSSAHash; hdp,smp: pointer; hlen,slen: word): boolean;
  {-Signature verification operation}
  { n,e   - signer's RSA public key pair (n, e) = (modulus, exponent)
    AHash - hash algorithm for signature
    hdp   - hash digest
    hlen  - length of hash digest
    smp   - signature message
    slen  - length of signature in bytes}
var
  k: word;
  s,m: mp_int;
  buf: pointer;
begin
  mp_pkcs1v15_verify := false;

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    {Check e here; n checked in mp_unsigned_bin_size}
    if mp_not_init(e) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pkcs1v15_verify');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if slen<>k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_verify: slen <> k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;

  {Get local memory for second encoded message}
  buf := mp_getmem(slen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_verify: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {Produce a second encoded message}
  if mp_pkcs1v15_emsa_encode(AHash, hdp, buf, hlen, slen) then begin
    mp_init2(s,m);
    if mp_error=MP_OKAY then begin
      {Convert the signature to an integer signature representative}
      mp_os2ip(smp, k, s);
      {Apply the RSA verification primitive}
      mp_rsavp(s,e,n,m);
      {Convert the second encode message to an integer representative}
      mp_os2ip(buf, k, s);
      {Verify that both integer representatives are the same}
      mp_pkcs1v15_verify := mp_is_eq(s,m);
      mp_clear2(s,m);
    end;
  end;
  fillchar(buf^,slen,0);
  mp_freemem(buf, slen);
end;


{---------------------------------------------------------------------------}
procedure mp_pkcs1v15_sign(const d, n: mp_int; AHash: TSSAHash; hdp,smp: pointer; hlen,smax: word; var slen: word);
  {-Sign a hash digest using RSA and EMSA-PKCS1-v1_5 encoding}
  { n,d   - signer's RSA private key pair (n, d) = (modulus, exponent)
    AHash - hash algorithm for signature
    hdp   - hash digest
    hlen  - length of hash digest
    smp   - signature message
    smax  - max length of signature buffer
    slen  - length of signature in bytes}
var
  k: word;
  m: mp_int;
  buf: pointer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    {Check d here; n checked in mp_unsigned_bin_size}
    if mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pkcs1v15_sign');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if smax<k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_sign: smax < k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  slen := k;
  {get local memory for encoded message}
  buf := mp_getmem(slen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_sign: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {initialize signature mp_int representative}
  mp_init(m);
  if mp_error=MP_OKAY then begin
    {Apply EMSA-PKCS1-v1_5 encoding: hash digest -> encoded message EM in buf}
    if not mp_pkcs1v15_emsa_encode(AHash, hdp, buf, hlen, k) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXUndef.Create('mp_pkcs1v15_sign: encode error');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$endif}
      set_mp_error(MP_UNDEF);
    end
    else begin
      {Convert the encoded message EM to mp_int representative m}
      mp_os2ip(buf, k, m);
      {Apply the RSA signature primitive to the message presentative m}
      mp_rsasp(m, d, n, m);
      {Convert the signature representative to octet string}
      mp_i2osp(m, smp, slen);
    end;
    mp_clear(m);
  end;
  fillchar(buf^,slen,0);
  mp_freemem(buf, slen);
end;


{---------------------------------------------------------------------------}
procedure mp_pkcs1v15_sign2(const prk: TPrivateKey; const n: mp_int; AHash: TSSAHash;
                            hdp,smp: pointer; hlen,smax: word; var slen: word);
  {-Sign a hash digest using RSA/CRT and EMSA-PKCS1-v1_5 encoding}
  { prk,n - signer's private key CRT record and modulus
    AHash - hash algorithm for signature
    hdp   - hash digest
    hlen  - length of hash digest
    smp   - signature message
    smax  - max length of signature buffer
    slen  - length of signature in bytes}
var
  k: word;
  m: mp_int;
  buf: pointer;
begin
  if mp_error<>MP_OKAY then exit;
  {Get size of modulus n in bytes}
  k := mp_unsigned_bin_size(n);
  if smax<k then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_pkcs1v15_sign2: smax < k');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
    exit;
  end;
  slen := k;
  {get local memory for encoded message}
  buf := mp_getmem(slen);
  if buf=nil then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMemory.Create('mp_pkcs1v15_sign2: buf=nil');
      {$else}
        RunError(MP_RTE_MEM);
      {$endif}
    {$else}
      set_mp_error(MP_MEM);
      exit;
    {$endif}
  end;

  {initialize signature mp_int representative}
  mp_init(m);
  if mp_error=MP_OKAY then begin
    {Apply EMSA-PKCS1-v1_5 encoding: hash digest -> encoded message EM in buf}
    if not mp_pkcs1v15_emsa_encode(AHash, hdp, buf, hlen, k) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXUndef.Create('mp_pkcs1v15_sign2: encode error');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$endif}
      set_mp_error(MP_UNDEF);
    end
    else begin
      {Convert the encoded message EM to mp_int representative m}
      mp_os2ip(buf, k, m);
      {Apply the RSA signature primitive to the message representative m}
      mp_rsasp2(m, prk, m);
      {Convert the signature representative to octet string}
      mp_i2osp(m, smp, slen);
    end;
    mp_clear(m);
  end;
  fillchar(buf^,slen,0);
  mp_freemem(buf, slen);
end;


end.

#! /bin/bash

# We need to diff gcc/libquadmath/math/ and glibc/sysdeps/ieee754/ldbl-128/

set gcc_dir=$HOME/gcc/libquadmath/math
set glibc_dir=$HOME/glibc/sysdeps/ieee754/ldbl-128

diff ${gcc_dir}/acoshq.c ${glibc_dir}/e_acoshl.c

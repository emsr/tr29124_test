#!  /bin/bash

tool="cp -f"

testdir="testsuite/special_functions"

rm -rf "testsuite"

mkdir -p $testdir/01_assoc_laguerre
mkdir -p $testdir/02_assoc_legendre
mkdir -p $testdir/03_beta
mkdir -p $testdir/04_comp_ellint_1
mkdir -p $testdir/05_comp_ellint_2
mkdir -p $testdir/06_comp_ellint_3
mkdir -p $testdir/07_cyl_bessel_i
mkdir -p $testdir/08_cyl_bessel_j
mkdir -p $testdir/09_cyl_bessel_k
mkdir -p $testdir/10_cyl_neumann
mkdir -p $testdir/11_ellint_1
mkdir -p $testdir/12_ellint_2
mkdir -p $testdir/13_ellint_3
mkdir -p $testdir/14_expint
mkdir -p $testdir/15_hermite
mkdir -p $testdir/16_laguerre
mkdir -p $testdir/17_legendre
mkdir -p $testdir/18_riemann_zeta
mkdir -p $testdir/19_sph_bessel
mkdir -p $testdir/20_sph_legendre
mkdir -p $testdir/21_sph_neumann

$tool testcase.h                $testdir

$tool check_assoc_laguerre.cc   $testdir/01_assoc_laguerre/check_value.cc
$tool check_assoc_legendre.cc   $testdir/02_assoc_legendre/check_value.cc
$tool check_beta.cc             $testdir/03_beta/check_value.cc
$tool check_comp_ellint_1.cc    $testdir/04_comp_ellint_1/check_value.cc
$tool check_comp_ellint_2.cc    $testdir/05_comp_ellint_2/check_value.cc
$tool check_comp_ellint_3.cc    $testdir/06_comp_ellint_3/check_value.cc
$tool check_cyl_bessel_i.cc     $testdir/07_cyl_bessel_i/check_value.cc
$tool check_cyl_bessel_j.cc     $testdir/08_cyl_bessel_j/check_value.cc
$tool check_cyl_bessel_k.cc     $testdir/09_cyl_bessel_k/check_value.cc
$tool check_cyl_neumann.cc      $testdir/10_cyl_neumann/check_value.cc
$tool check_ellint_1.cc         $testdir/11_ellint_1/check_value.cc
$tool check_ellint_2.cc         $testdir/12_ellint_2/check_value.cc
$tool check_ellint_3.cc         $testdir/13_ellint_3/check_value.cc
$tool check_expint_neg.cc       $testdir/14_expint/check_value_neg.cc
$tool check_expint_pos.cc       $testdir/14_expint/check_value_pos.cc

$tool check_laguerre.cc         $testdir/16_laguerre/check_value.cc
$tool check_legendre.cc         $testdir/17_legendre/check_value.cc
$tool check_riemann_zeta_neg.cc $testdir/18_riemann_zeta/check_value_neg.cc
$tool check_riemann_zeta_pos.cc $testdir/18_riemann_zeta/check_value_pos.cc
$tool check_sph_bessel.cc       $testdir/19_sph_bessel/check_value.cc
$tool check_sph_legendre.cc     $testdir/20_sph_legendre/check_value.cc
$tool check_sph_neumann.cc      $testdir/21_sph_neumann/check_value.cc

{T_SFX_M - regression test main unit for  SPECFUNX  (c) 2010-2018  W.Ehrhardt}

unit t_sfx_m;

interface

{$i STD.INC}


{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


procedure test_sfx_main;
  {-SpecFunX regression test via single procedure call}


implementation

uses
  {$ifdef debug}
    sfBasic,
  {$endif}
  specfunx, t_sfx0,
  t_sfx1,t_sfx1a,t_sfx1b,
  t_sfx2,t_sfx2a,t_sfx2b,t_sfx2c,
  t_sfx3,t_sfx3a,t_sfx3b,t_sfx3c,t_sfx3d,
  t_sfx4,t_sfx4a,
  t_sfx5,t_sfx5a,
  t_sfx6,t_sfx6a,t_sfx6b,
  t_sfx7,t_sfx7a,
  t_sfx8,t_sfx8a,t_sfx8b,
  t_sfx9,t_sfx9a;


{---------------------------------------------------------------------------}
procedure test_sfx_main;
  {-SpecFunX regression test via single procedure call}

label
  done;

begin
  {$ifdef debug}
    sfc_debug_output := true;
  {$endif}
  writeln('Test SpecFunX V', SpecFunX_Version, ' [extended precision]  -  (c) 2010-2018 W.Ehrhardt');

  total_cnt    := 0;
  total_failed := 0;

{
  test_kelvin_derx;
  goto done;
}
  test_bessel_i0x;
  test_bessel_i0ex;
  test_bessel_i1x;
  test_bessel_i1ex;
  test_bessel_j0x;
  test_bessel_j1x;
  test_bessel_k0x;
  test_bessel_k0ex;
  test_bessel_k1x;
  test_bessel_k1ex;
  test_bessel_y0x;
  test_bessel_y1x;

  test_bessel_jvx;
  test_bessel_yvx;
  test_bessel_ivx;
  test_bessel_kvx;
  test_bessel_kvex;
  test_bessel_ivex;
  test_jnynx;
  test_inknx;
  test_bessel_lambdax;

  test_j0_intx;
  test_y0_intx;
  test_i0_intx;
  test_k0_intx;

  test_sph_bess_jnx;
  test_sph_bess_ynx;
  test_sph_bessel_inx;
  test_sph_bessel_inex;
  test_sph_bessel_knx;
  test_sph_bessel_knex;

  test_airy_aix;
  test_airy_bix;
  test_airy_aipx;
  test_airy_bipx;
  test_airy_aisx;
  test_airy_bisx;
  test_airy_scorerx;

  test_berx;
  test_beix;
  test_kerx;
  test_keix;
  test_kelvin_derx;

  test_struve_h0x;
  test_struve_h1x;
  test_struve_l0x;
  test_struve_l1x;
  test_struve_hx;
  test_struve_lx;

  test_CoulombCLx;
  test_CoulombSLx;
  test_CoulombFFpx;
  test_CoulombGGpx;

  test_synchfx;
  test_synchgx;

  {t_sfx1}
  test_facx;
  test_gammax;
  test_igammax;
  test_igammalx;
  test_igammanx;
  test_igammatx;
  test_gammastarx;
  test_gamma1pm1x;
  test_gamma_ratiox;
  test_pochhammerx;
  test_poch1x;
  test_binomialx;
  test_lnbinomialx;
  test_dfacx;
  test_lnfacx;
  test_rgammax;
  test_inv_gammax;
  test_igammap_derx;

  test_zetax;
  test_zeta1px;
  test_zetaintx;
  test_igamma_invx;
  test_incgammax;
  test_zetam1x;
  test_primezetax;
  test_etam1x;
  test_etax;
  test_polylogx;
  test_polylogrx;
  test_bose_einsteinx;
  test_fermi_diracx;
  test_fermi_dirac_rx;
  test_fermi_dirac_halfx;
  test_LerchPhix;
  test_DirichletBetax;
  test_DirichletLambdax;
  test_LegendreChix;

  test_zetahx;
  test_harmonicx;
  test_harmonic2x;
  test_lngammax;
  test_lngamma_invx;
  test_lngamma1px;
  test_psix;
  test_psistarx;
  test_psi_invx;
  test_trigammax;
  test_tetrapentagammax;
  test_polygammax;
  test_lnBarnesGx;
  test_BatemanGx;
  test_lnbetax;
  test_betax;
  test_beta3x;
  test_ibetax;
  test_ibeta_invx;
  test_lobachevsky_cx;
  test_lobachevsky_sx;

  test_carlsonx;
  test_bulirschx;
  test_maplex;
  test_mathematicax;
  test_legendrex;
  test_comp_ellintx;
  test_EllipticKimx;
  test_EllipticECimx;
  test_EllipticPiCimx;
  test_heuman_lambdax;
  test_jacobi_zetax;

  test_jacobiPQx;
  test_jacobi_amx;
  test_sncndnx;
  test_theta1px;
  test_theta2x;
  test_theta3x;
  test_theta4x;
  test_jacobi_thetax;
  test_jacobi_thetax_relations;
  test_nevillex;
  test_EllipticModulusx;
  test_EllipticNomex;
  test_lemniscatex;

  test_jacobi_arcsnx;
  test_jacobi_arccnx;
  test_jacobi_arcdnx;
  test_jacobi_arcscx;
  test_jacobi_arccsx;
  test_jacobi_arcncx;
  test_jacobi_arcnsx;
  test_jacobi_arcndx;
  test_jacobi_arccdx;
  test_jacobi_arcdcx;
  test_jacobi_arcsdx;
  test_jacobi_arcdsx;

  test_detaix;
  test_emlambdax;
  test_KleinJx;
  test_wplx;
  test_wpex;
  test_wpe_derx;
  test_wpe_invx;
  test_wpgx;
  test_wpg_derx;
  test_wpg_invx;

  test_agmx;
  test_cl2x;
  test_dilogx;
  test_trilogx;
  test_ti2x;
  test_tix;
  test_lambertwx;
  test_lambertw1x;
  test_debyex;
  test_einsteinx;
  test_RiemannRx;
  test_RiemannR_invx;
  test_bernpolyx;
  test_besselpolyx;
  test_cosintx;
  test_sinintx;
  test_fibfunx;
  test_fibpolyx;
  test_lucpolyx;
  test_catalanx;
  test_eulerx;
  test_eulerpolyx;
  test_euler_qx;
  test_keplerx;
  test_langevinx;
  test_langevin_invx;
  test_bringx;
  test_omegax;
  test_transportx;
  test_rrcfx;
  test_exprelnx;
  test_expnx;


  test_erfx;
  test_erfcx;
  test_erfcex;
  test_erfgx;
  test_inerfcx;
  test_erfix;
  test_erfh2x;
  test_dawsonx;
  test_erfinvx;
  test_dawson2x;
  test_erfcinvx;
  test_erfce_invx;
  test_erfi_invx;
  test_expint3x;
  test_fresnelx;
  test_fresnelfgx;
  test_gsix;
  test_erf_zx;
  test_erf_px;
  test_erf_qx;
  test_marcumqx;
  test_owentx;

  test_beta_cdfx;
  test_beta_invx;
  test_beta_pdfx;
  test_cauchy_cdfx;
  test_cauchy_pdfx;
  test_cauchy_invx;
  test_evt1_cdfx;
  test_evt1_pdfx;
  test_evt1_invx;
  test_exp_pdfx;
  test_exp_cdfx;
  test_exp_invx;
  test_f_pdfx;
  test_f_cdfx;
  test_f_invx;
  test_laplace_cdfx;
  test_laplace_invx;
  test_laplace_pdfx;
  test_levy_cdfx;
  test_levy_invx;
  test_levy_pdfx;
  test_logistic_cdfx;
  test_logistic_invx;
  test_logistic_pdfx;
  test_lognormal_cdfx;
  test_lognormal_invx;
  test_lognormal_pdfx;
  test_negbinom_pmfx;
  test_negbinom_cdfx;
  test_normal_pdfx;
  test_normal_cdfx;
  test_normal_invx;
  test_normstd_pdfx;
  test_normstd_cdfx;
  test_normstd_invx;
  test_pareto_cdfx;
  test_pareto_invx;
  test_pareto_pdfx;
  test_triangular_pdfx;
  test_triangular_invx;
  test_triangular_cdfx;
  test_uniform_cdfx;
  test_uniform_invx;
  test_uniform_pdfx;
  test_weibull_cdfx;
  test_weibull_invx;
  test_weibull_pdfx;
  test_binomial_pmfx;
  test_binomial_cdfx;
  test_poisson_pmfx;
  test_poisson_cdfx;
  test_hypergeo_pmfx;
  test_hypergeo_cdfx;
  test_rayleigh_pdfx;
  test_rayleigh_cdfx;
  test_rayleigh_invx;
  test_maxwell_cdfx;
  test_maxwell_pdfx;
  test_moyal_pdfx;
  test_moyal_cdfx;
  test_moyal_invx;
  test_kolmogorov_cdfx;
  test_kolmogorov_invx;
  test_kumaraswamy_pdfx;
  test_kumaraswamy_cdfx;
  test_kumaraswamy_invx;
  test_nakagami_cdfx;
  test_nakagami_pdfx;
  test_nakagami_invx;
  test_zipf_pmfx;
  test_zipf_cdfx;
  test_logseries_cdfx;
  test_logseries_pmfx;
  test_wald_cdfx;
  test_wald_invx;
  test_wald_pdfx;

  test_chi2_pdfx;
  test_chi2_cdfx;
  test_chi2_invx;
  test_chi_distx;
  test_gamma_invx;
  test_gamma_pdfx;
  test_gamma_cdfx;
  test_invgamma_pdfx;
  test_invgamma_cdfx;
  test_invgamma_invx;
  test_t_cdfx;
  test_t_pdfx;
  test_t_invx;
  test_maxwell_invx;

  test_e1x;
  test_e1sx;
  test_eix;
  test_einx;
  test_eisx;
  test_eisx2x;
  test_ei_invx;
  test_lix;
  test_li_invx;
  test_enx;
  test_eibetax;
  test_geix;
  test_chix;
  test_shix;
  test_cix;
  test_six;
  test_ssix;
  test_cinx;
  test_cinhx;

  test_chebyshev_tx;
  test_chebyshev_ux;
  test_gegenbauer_cx;
  test_chebyshev_vx;
  test_chebyshev_wx;
  test_chebfunx;
  test_hermite_hx;
  test_hermite_hex;
  test_jacobi_px;
  test_laguerrex;
  test_legendre_plx;
  test_legendre_qlx;
  test_legendre_plmx;
  test_legendre_qlmx;
  test_toroidal_qlmx;
  test_toroidal_plmx;
  test_zernike_rx;
  test_spherical_harmonix;
  test_orthopolyx;

  test_hyperg_1F1x;
  test_hyperg_1F1rx;
  test_hyperg_ux;
  test_hyperg_2F1x;
  test_hyperg_2F1rx;
  test_whittakerx;
  test_hyperg_0F1x;
  test_hyperg_0F1rx;
  test_hyperg_2F0x;
  test_cylinderdx;
  test_cylinderux;
  test_cylindervx;
  test_hermitehx;


done:

  if total_failed>0 then begin
    writeln('***** Total number of failed tests: ',total_failed,' of ',total_cnt);
  end
  else writeln('Passed. All ',total_cnt, ' test were OK.');
end;

end.


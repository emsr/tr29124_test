{T_SFD_M - regression test main unit for SPECFUN  (c) 2010-2015  W.Ehrhardt}

unit t_sfd_m;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

{$ifdef NOBASM}
{$undef BASM}
{$endif}

procedure test_sfd_main;
  {-SpecFun regression test via single procedure call}

implementation

uses
  {$ifdef debug}
    sfBasic,
  {$endif}
  specfun, t_sfd0,
  t_sfd1,t_sfd1a,t_sfd1b,
  t_sfd2,t_sfd2a,t_sfd2b,
  t_sfd3,t_sfd3a,t_sfd3b,t_sfd3c,t_sfd3d,
  t_sfd4,t_sfd4a,t_sfd5,t_sfd5a,t_sfd6,t_sfd6a,
  t_sfd7,t_sfd8,t_sfd8a,t_sfd9,t_sfd9a;


{---------------------------------------------------------------------------}
procedure test_sfd_main;
  {-SpecFun regression test via single procedure call}
label
  done;

begin
  {$ifdef debug}
    sfc_debug_output := true;
  {$endif}
  writeln('Test SpecFun V', SpecFun_Version, ' [double precision]  -  (c) 2010-2015 W.Ehrhardt');

  total_cnt    := 0;
  total_failed := 0;

{
  test_lemniscate;
  test_cylinderv;
  test_polygamma;
  goto done;
}
  test_bessel_i0;
  test_bessel_i0e;
  test_bessel_i1;
  test_bessel_i1e;
  test_bessel_j0;
  test_bessel_j1;
  test_bessel_k0;
  test_bessel_k0e;
  test_bessel_k1;
  test_bessel_k1e;
  test_bessel_y0;
  test_bessel_y1;

  test_bessel_jv;
  test_bessel_yv;
  test_bessel_iv;
  test_bessel_kv;
  test_bessel_ive;
  test_bessel_kve;
  test_jnyn;
  test_inkn;

  test_sph_bess_jn;
  test_sph_bess_yn;
  test_sph_bessel_in;
  test_sph_bessel_ine;
  test_sph_bessel_kn;
  test_sph_bessel_kne;

  test_airy_ai;
  test_airy_bi;
  test_airy_aip;
  test_airy_bip;
  test_airy_ais;
  test_airy_bis;
  test_airy_scorer;

  test_ber;
  test_bei;
  test_ker;
  test_kei;

  test_struve_h0;
  test_struve_h1;
  test_struve_l0;
  test_struve_l1;
  test_struve_h;
  test_struve_l;

  test_gamma;
  test_gamma1pm1;
  test_gammastar;
  test_gamma_ratio;
  test_pochhammer;
  test_poch1;
  test_binomial;
  test_fac;
  test_dfac;
  test_lnfac;
  test_lngamma;
  test_lngamma1p;
  test_lngamma_inv;
  test_rgamma;
  test_incgamma;
  test_igamma;
  test_igammal;
  test_igammat;
  test_beta3;

  test_zeta;
  test_zeta1p;
  test_zetaint;
{$ifndef VER50}
  test_igamma_inv;
  test_zetam1;
  test_primezeta;
  test_etam1;
  test_eta;
  test_polylog;
  test_polylogr;
  test_fermi_dirac;
  test_fermi_dirac_half;
  test_LerchPhi;
  test_DirichletBeta;
  test_DirichletLambda;
  test_LegendreChi;
{$endif}
  test_zetah;
  test_harmonic;
  test_harmonic2;
  test_psi;
  test_psi_inv;
  test_trigamma;
  test_tetrapentagamma;
  test_polygamma;
  test_BatemanG;
  test_beta;
  test_ibeta;
  test_ibeta_inv;
  test_lnbeta;
  test_lobachevsky_c;
  test_lobachevsky_s;

{$ifndef VER50}
  test_carlson;
{$endif}
  test_bulirsch;
  test_maple;
  test_legendre;
  test_comp_ellint;
  test_EllipticKim;
  test_EllipticECim;
  test_heuman_lambda;
  test_jacobi_zeta;

  test_jacobiPQ;
  test_jacobi_am;
  test_sncndn;
  test_theta1p;
  test_theta2;
  test_theta3;
  test_theta4;
  test_jacobi_theta;
  test_jacobi_theta_relations;
  test_EllipticModulus;
  test_EllipticNome;
  test_lemniscate;

  test_jacobi_arcsn;
  test_jacobi_arccn;
  test_jacobi_arcdn;
  test_jacobi_arcsc;
  test_jacobi_arccs;
  test_jacobi_arcnc;
  test_jacobi_arcns;
  test_jacobi_arcnd;
  test_jacobi_arccd;
  test_jacobi_arcdc;
  test_jacobi_arcsd;
  test_jacobi_arcds;

  test_agm;
  test_cl2;
  test_dilog;
  test_trilog;
  test_sncndn;
  test_ti2;
  test_lambertw;
  test_lambertw1;
  test_li_inv;
  test_debye;
  test_RiemannR;
  test_bernpoly;
  test_cosint;
  test_sinint;
  test_fibpoly;
  test_lucpoly;
  test_catalan;
  test_euler;
  test_kepler;

  test_erf;
  test_erfc;
  test_erfce;
  test_erfg;
  test_inerfc;
  test_dawson;
  test_dawson2;
  test_erfi;
  test_erfinv;
  test_erfcinv;
  test_expint3;
  test_fresnel;
  test_gsi;
  test_erf_z;
  test_erf_p;
  test_erf_q;

  test_beta_cdf;
  test_beta_inv;
  test_beta_pdf;
  test_cauchy_cdf;
  test_cauchy_inv;
  test_cauchy_pdf;
  test_evt1_cdf;
  test_evt1_pdf;
  test_evt1_inv;
  test_exp_cdf;
  test_exp_inv;
  test_exp_pdf;
  test_f_pdf;
  test_f_cdf;
  test_f_inv;
  test_laplace_cdf;
  test_laplace_inv;
  test_laplace_pdf;
  test_levy_cdf;
  test_levy_inv;
  test_levy_pdf;
  test_logistic_cdf;
  test_logistic_inv;
  test_logistic_pdf;
  test_lognormal_cdf;
  test_lognormal_inv;
  test_lognormal_pdf;
  test_negbinom_pmf;
  test_negbinom_cdf;
  test_normal_pdf;
  test_normal_cdf;
  test_normal_inv;
  test_normstd_pdf;
  test_normstd_cdf;
  test_normstd_inv;
  test_pareto_cdf;
  test_pareto_inv;
  test_pareto_pdf;
  test_triangular_pdf;
  test_triangular_inv;
  test_triangular_cdf;
  test_uniform_cdf;
  test_uniform_inv;
  test_uniform_pdf;
  test_weibull_cdf;
  test_weibull_inv;
  test_weibull_pdf;
  test_binomial_pmf;
  test_binomial_cdf;
  test_poisson_pmf;
  test_poisson_cdf;
  test_hypergeo_pmf;
  test_hypergeo_cdf;
  test_rayleigh_pdf;
  test_rayleigh_cdf;
  test_rayleigh_inv;
  test_maxwell_cdf;
  test_maxwell_pdf;
  test_kumaraswamy_pdf;
  test_kumaraswamy_cdf;
  test_kumaraswamy_inv;
  test_moyal_pdf;
  test_moyal_cdf;
  test_moyal_inv;
  test_zipf_pmf;
  test_zipf_cdf;
  test_logseries_cdf;
  test_logseries_pmf;
  test_wald_cdf;
  test_wald_inv;
  test_wald_pdf;

{$ifndef VER50}
  test_chi2_pdf;
  test_chi2_cdf;
  test_chi2_inv;
  test_gamma_inv;
  test_gamma_pdf;
  test_gamma_cdf;
  test_invgamma_pdf;
  test_invgamma_cdf;
  test_invgamma_inv;
  test_t_cdf;
  test_t_pdf;
  test_t_inv;
  test_maxwell_inv;
{$endif}

  test_e1;
  test_ei;
  test_gei;
  test_ein;
  test_ei_inv;
  test_li;
  test_en;
  test_chi;
  test_shi;
  test_ci;
  test_cin;
  test_cinh;
  test_si;
  test_ssi;

  test_chebyshev_t;
  test_chebyshev_u;
  test_chebyshev_v;
  test_chebyshev_w;
  test_gegenbauer_c;
  test_hermite_h;
  test_jacobi_p;
  test_laguerre;
  test_legendre_pl;
  test_legendre_ql;
  test_legendre_plm;
  test_legendre_qlm;
  test_toroidal_plm;
  test_toroidal_qlm;
  test_zernike_r;
  test_spherical_harmonics;
  test_orthopoly;

{$ifndef VER50}
  test_hyperg_1F1;
  test_hyperg_1F1r;
  test_hyperg_u;
  test_hyperg_2F1;
  test_hyperg_2F1r;
  test_whittaker;
  test_hyperg_0F1;
  test_hyperg_0F1r;
  test_cylinderd;
  test_cylinderu;
  test_cylinderv;
  test_hermiteh;
{$endif}


done:

  if total_failed>0 then begin
    writeln('***** Total number of failed tests: ',total_failed,' of ',total_cnt);
  end
  else writeln('Passed. All ',total_cnt, ' test were OK.');
end;


end.

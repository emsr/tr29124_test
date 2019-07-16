{Part 3c of regression test for SPECFUN unit  (c) 2011  W.Ehrhardt}

unit t_sfd3c;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_binomial_pmf;
procedure test_binomial_cdf;

procedure test_negbinom_pmf;
procedure test_negbinom_cdf;

procedure test_poisson_pmf;
procedure test_poisson_cdf;

procedure test_hypergeo_pmf;
procedure test_hypergeo_cdf;

procedure test_rayleigh_cdf;
procedure test_rayleigh_pdf;
procedure test_rayleigh_inv;

procedure test_maxwell_cdf;
procedure test_maxwell_pdf;
procedure test_maxwell_inv;

procedure test_evt1_cdf;
procedure test_evt1_pdf;
procedure test_evt1_inv;

procedure test_kumaraswamy_pdf;
procedure test_kumaraswamy_cdf;
procedure test_kumaraswamy_inv;

procedure test_kolmogorov_cdf;
procedure test_kolmogorov_inv;

implementation

uses
  specfun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_binomial_pmf;
var
  y,f,p: double;
  n: longint;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','binomial_pmf');

  {statevalf[pf,binomiald[n,p]](k);}

  p := 0.75;
  n := 40;

  y := binomial_pmf(p,n,-1);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,0);
  f := 0.8271806125530276749e-24;
  testrel(2, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,1);
  f := 0.9926167350636332098e-22;
  testrel(3, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,5);
  f := 0.1322628249026646101e-15;
  testrel(4, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,10);
  f := 0.4140329018188032423e-10;
  testrel(5, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,20);
  f := 0.3975770213715739165e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,30);
  f := 0.1443643463562567674;
  testrel(7, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,35);
  f := 0.2723174275324595051e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,39);
  f := 0.1340878021551666236e-3;
  testrel(9, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,40);
  f := 0.1005658516163749677e-4;
  testrel(10, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,41);
  f := 0.0;
  testrel(11, NE, y, f, cnt,failed);


  p := 0.03125;
  n := 1234;
  y := binomial_pmf(p,n,-1);
  f := 0.0;
  testrel(12, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,0);
  f := 0.9666225996407504662e-17;
  testrel(13, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,1);
  f := 0.3847781574053826049e-15;
  testrel(14, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,10);
  f := 0.2565372608556586559e-7;
  testrel(15, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,20);
  f := 0.3396080149874057448e-3;
  testrel(16, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,40);
  f := 0.6226731743699481439e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,60);
  f := 0.2685931148682441038e-3;
  testrel(18, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,200);
  f := 0.4615245266363444125e-79;
  testrel(19, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,1000);
  f := 0.2353735912995348529e-1249;
  testrel(20, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,n-1);
  f := 0.1688898720849300065e-1852;
  testrel(21, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,n);
  f := 0.4414959797274272141e-1857;
  testrel(22, NE, y, f, cnt,failed);

  y := binomial_pmf(p,n,n+1);
  f := 0.0;
  testrel(23, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_binomial_cdf;
var
  y,f,p: double;
  n: longint;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','binomial_cdf');

  {statevalf[dcdf,binomiald[n,p]](k);}
  p := 0.75;
  n := 40;

  y := binomial_cdf(p,n,-1);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,0);
  f := 0.8271806125530276749e-24;
  testrel(2, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,1);
  f := 0.1000888541189163487e-21;
  testrel(3, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,5);
  f := 0.1386126694303024141e-15;
  testrel(4, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,10);
  f := 0.4630880897739289242e-10;
  testrel(5, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,20);
  f := 0.5724311071761385043e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,30);
  f := 0.5604602683274669684;
  testrel(7, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,35);
  f := 0.9839577601812303383;
  testrel(8, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,39);
  f := 0.9999899434148383625;
  testrel(9, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,40);
  f := 1.0;
  testrel(10, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,41);
  f := 1.0;
  testrel(11, NE, y, f, cnt,failed);


  p := 0.03125;
  n := 1234;

  y := binomial_cdf(p,n,-1);
  f := 0.0;
  testrel(12, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,0);
  f := 0.9666225996407504662e-17;
  testrel(13, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,1);
  f := 0.3944443834017901096e-15;
  testrel(14, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,10);
  f := 0.3398353413634645458e-7;
  testrel(15, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,20);
  f := 0.6626922527501145631e-3;
  testrel(16, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,40);
  f := 0.6330323277786663291;
  testrel(17, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,60);
  f := 0.9995867214272379805;
  testrel(18, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,70);
  f := 0.9999988156877154370;
  testrel(19, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,80);
  f := 0.9999999991849552041;
  testrel(20, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,100);
  f := 0.9999999999999999907;
  testrel(21, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,1000);
  f := 1.0;
  testrel(22, NE, y, f, cnt,failed);

  y := binomial_cdf(p,n,1001);
  f := 1.0;
  testrel(23, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_poisson_pmf;
var
  y,f,mu: double;
  k: longint;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','poisson_pmf');

  {statevalf[pf,poisson[mu]](k);}

  mu := 7.5;

  k := 0;
  f := 0.5530843701478335831e-3;
  y := poisson_pmf(mu,k);
  testrel(1, NE, y, f, cnt,failed);

  k := 1;
  f := 0.4148132776108751873e-2;
  y := poisson_pmf(mu,k);
  testrel(2, NE, y, f, cnt,failed);

  k := 2;
  f := 0.1555549791040781952e-1;
  y := poisson_pmf(mu,k);
  testrel(3, NE, y, f, cnt,failed);

  k := 5;
  f := 0.1093745946825549810;
  y := poisson_pmf(mu,k);
  testrel(4, NE, y, f, cnt,failed);

  k := 8;
  f := 0.1373285926538776269;
  y := poisson_pmf(mu,k);
  testrel(5, NE, y, f, cnt,failed);

  k := 10;
  f := 0.8583037040867351678e-1;
  y := poisson_pmf(mu,k);
  testrel(6, NE, y, f, cnt,failed);

  k := 20;
  f := 0.7209282379462169336e-4;
  y := poisson_pmf(mu,k);
  testrel(7, NE, y, f, cnt,failed);

  k := 40;
  f := 0.6817055868620604552e-16;
  y := poisson_pmf(mu,k);
  testrel(8, NE, y, f, cnt,failed);

  k := 50;
  f := 0.1029863539144823660e-23;
  y := poisson_pmf(mu,k);
  testrel(9, NE, y, f, cnt,failed);

  mu := 1000.0;

  k := 0;
  f := 0.5075958897549456765e-434;
  y := poisson_pmf(mu,k);
  testrel(10, NE, y, f, cnt,failed);
  k := 100;
  f := 0.5438942180826245858e-292;
  y := poisson_pmf(mu,k);
  testrel(11, NE, y, f, cnt,failed);

  k := 800;
  f := 0.6583151641880508578e-11;
  y := poisson_pmf(mu,k);
  testrel(12, NE, y, f, cnt,failed);

  k := 900;
  f := 0.7516954352125952229e-4;
  y := poisson_pmf(mu,k);
  testrel(13, NE, y, f, cnt,failed);

  k := 1000;
  f := 0.1261461134872149972e-1;
  y := poisson_pmf(mu,k);
  testrel(14, NE, y, f, cnt,failed);

  k := 1100;
  f := 0.9498944242299507576e-4;
  y := poisson_pmf(mu,k);
  testrel(15, NE, y, f, cnt,failed);

  k := 1200;
  f := 0.7992642848843570799e-10;
  y := poisson_pmf(mu,k);
  testrel(16, NE, y, f, cnt,failed);

  k := 1300;
  f := 0.1606560638609706058e-19;
  y := poisson_pmf(mu,k);
  testrel(17, NE, y, f, cnt,failed);

  k := 1400;
  f := 0.1466773344196568303e-32;
  y := poisson_pmf(mu,k);
  testrel(18, NE, y, f, cnt,failed);

  y := poisson_pmf(1000000,1000000-5000);
  f := 0.14596440994146676390e-8;
  testrel(19, NE, y, f, cnt,failed);

  y := poisson_pmf(1000000,1000000-1);
  f := 0.39894224715624402976e-3;
  testrel(20, NE, y, f, cnt,failed);

  y := poisson_pmf(1000000,1000000);
  f := 0.39894224715624402970e-3;
  testrel(21, NE, y, f, cnt,failed);

  y := poisson_pmf(1000000,1000000+1);
  f := 0.3989418482143958153e-3;
  testrel(22, NE, y, f, cnt,failed);

  y := poisson_pmf(1000000,1000000+5000);
  f := 0.1514158102861422074e-8;
  testrel(23, NE, y, f, cnt,failed);

  y := poisson_pmf(1e-10,2);
  f := 0.4999999999500000000025e-20;
  testrel(24, NE, y, f, cnt,failed);

  y := poisson_pmf(1e-19,2);
  f := 0.49999999999999999995e-38;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_poisson_cdf;
var
  y,f,mu: double;
  k: longint;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','poisson_cdf');

  {statevalf[dcdf,poisson[mu]](k);}

  mu := 7.5;

  k := 0;
  f := 0.5530843701478335831e-3;
  y := poisson_cdf(mu,k);
  testrel(1, NE, y, f, cnt,failed);

  k := 1;
  f := 0.4701217146256585456e-2;
  y := poisson_cdf(mu,k);
  testrel(2, NE, y, f, cnt,failed);

  k := 2;
  f := 0.2025671505666440498e-1;
  y := poisson_cdf(mu,k);
  testrel(3, NE, y, f, cnt,failed);

  k := 5;
  f := 0.2414364509702755889;
  y := poisson_cdf(mu,k);
  testrel(4, NE, y, f, cnt,failed);

  k := 8;
  f := 0.6619671191414830773;
  y := poisson_cdf(mu,k);
  testrel(5, NE, y, f, cnt,failed);

  k := 10;
  f := 0.8622379834283879498;
  y := poisson_cdf(mu,k);
  testrel(6, NE, y, f, cnt,failed);

  k := 20;
  f := 0.999961343358912791;
  y := poisson_cdf(mu,k);
  testrel(7, NE, y, f, cnt,failed);

  k := 40;
  f := 0.999999999999999985;
  y := poisson_cdf(mu,k);
  testrel(8, NE, y, f, cnt,failed);

  k := 50;
  f := 1.0;
  y := poisson_cdf(mu,k);
  testrel(9, NE, y, f, cnt,failed);

  mu := 1000.0;

  k := 0;
  f := 0.5075958897549456765e-434;
  y := poisson_cdf(mu,k);
  testrel(10, NE, y, f, cnt,failed);

  k := 100;
  f := 0.6042524933789373681e-292;
  y := poisson_cdf(mu,k);
  testrel(11, NE, y, f, cnt,failed);
  k := 800;
  f := 0.3229888722729021536e-10;
  y := poisson_cdf(mu,k);
  testrel(12, NE, y, f, cnt,failed);

  k := 900;
  f := 0.6977673277963067821e-3;
  y := poisson_cdf(mu,k);
  testrel(13, NE, y, f, cnt,failed);

  k := 1000;
  f := 0.5084093671685059912;
  y := poisson_cdf(mu,k);
  testrel(14, NE, y, f, cnt,failed);

  k := 1100;
  f := 0.9991323590365564379;
  y := poisson_cdf(mu,k);
  testrel(15, NE, y, f, cnt,failed);

  k := 1200;
  f := 0.9999999996115060429;
  y := poisson_cdf(mu,k);
  testrel(16, NE, y, f, cnt,failed);

  k := 1300;
  f := 0.9999999999999999999;
  y := poisson_cdf(mu,k);
  testrel(17, NE, y, f, cnt,failed);

  k := 1400;
  f := 1.0;
  y := poisson_cdf(mu,k);
  testrel(18, NE, y, f, cnt,failed);

  y := poisson_cdf(1000000,1000000-5000);
  f := 0.2814820383896531442e-6;
  testrel(19, NE, y, f, cnt,failed);

  y := poisson_cdf(1000000,1000000-1);
  f := 0.4998670192391274088;
  testrel(20, NE, y, f, cnt,failed);

  y := poisson_cdf(1000000,1000000);
  f := 0.5002659614862836528;
  testrel(21, NE, y, f, cnt,failed);

  y := poisson_cdf(1000000,1000000+1);
  f := 0.5006649033344980486;
  testrel(22, NE, y, f, cnt,failed);

  y := poisson_cdf(1000000,1000000+5000);
  f := 0.9999997081107532997;
  testrel(23, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_negbinom_pmf;
var
  y,f,p,r: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','negbinom_pmf');

  {statevalf[pf,negativebinomial[r,p]](k);}

  p := 0.25;
  r := 12.5;

  y := negbinom_pmf(p,r,-1);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,0);
  f := 0.298023223876953125e-7;
  testrel(2, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,5);
  f := 0.3688099297960434342e-4;
  testrel(3, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,10);
  f := 0.8053113962378047894e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,20);
  f := 0.1315307743701173527e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,30);
  f := 0.3169313833290691748e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,40);
  f := 0.2995260427556122389e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,50);
  f := 0.1617482549671340422e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,100);
  f := 0.1394806553692310162e-4;
  testrel(9, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,1000);
  f := 0.8518841783313850163e-106;
  testrel(10, NE, y, f, cnt,failed);

  {case Pascal distribution: integer r}
  p := 0.5;
  r := 70;

  y := negbinom_pmf(p,r,-1);
  f := 0.0;
  testrel(11, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,0);
  f := 0.8470329472543003391e-21;
  testrel(12, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,15);
  f := 0.3828491912209976695e-9;
  testrel(13, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,30);
  f := 0.1621948340609462764e-4;
  testrel(14, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,60);
  f := 0.2566661064111152397e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,70);
  f := 0.3365662227431640229e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,80);
  f := 0.2179054458117639441e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,100);
  f := 0.1786171519214098839e-2;
  testrel(18, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,200);
  f := 0.9637265809589711395e-16;
  testrel(19, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,500);
  f := 0.2601453029343685356e-81;
  testrel(20, NE, y, f, cnt,failed);


  {case geometric distribution: r=1}
  p := 0.75;
  r := 1;

  y := negbinom_pmf(p,r,-1);
  f := 0.0;
  testrel(21, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,0);
  f := 0.75;
  testrel(22, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,1);
  f := 0.1875;
  testrel(23, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,2);
  f := 0.46875e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,4);
  f := 0.29296875e-2;
  testrel(25, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,6);
  f := 0.18310546875e-3;
  testrel(26, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,8);
  f := 0.11444091796875e-4;
  testrel(27, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,10);
  f := 0.7152557373046875e-6;
  testrel(28, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,25);
  f := 0.6661338147750939243e-15;
  testrel(29, NE, y, f, cnt,failed);

  y := negbinom_pmf(p,r,50);
  f := 0.5916456789157588541e-30;
  testrel(30, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_negbinom_cdf;
var
  y,f,p,r: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','negbinom_cdf');

  {statevalf[pf,negativebinomial[r,p]](k);}

  p := 0.25;
  r := 12.5;

  y := negbinom_cdf(p,r,-1);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,0);
  f := 0.298023223876953125e-7;
  testrel(2, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,5);
  f := 0.5863341687017964432e-4;
  testrel(3, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,10);
  f := 0.1896295783055346033e-2;
  testrel(4, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,20);
  f := 0.6270215545434477773e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,30);
  f := 0.3053534137005089327;
  testrel(6, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,40);
  f := 0.6309566199205073919;
  testrel(7, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,50);
  f := 0.8561377403191092259;
  testrel(8, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,100);
  f := 0.9999310806969732651;
  testrel(9, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,250);
  f := 0.9999999999999999998;
  testrel(10, NE, y, f, cnt,failed);

  {case Pascal distribution: integer r}
  p := 0.5;
  r := 70;

  y := negbinom_cdf(p,r,-1);
  f := 0.0;
  testrel(11, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,0);
  f := 0.8470329472543003391e-21;
  testrel(12, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,15);
  f := 0.5862470113599904290e-9;
  testrel(13, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,30);
  f := 0.3925069822796834811e-4;
  testrel(14, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,60);
  f := 0.2150104419876710213;
  testrel(15, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,70);
  f := 0.5336566222743164023;
  testrel(16, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,80);
  f := 0.8154192651691056915;
  testrel(17, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,90);
  f := 0.9517153233041100212;
  testrel(18, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,100);
  f := 0.9914156513274205386;
  testrel(19, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,200);
  f := 0.9999999999999998044;
  testrel(20, NE, y, f, cnt,failed);


  {case geometric distribution: r=1}
  p := 0.75;
  r := 1;

  y := negbinom_cdf(p,r,-1);
  f := 0.0;
  testrel(21, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,0);
  f := 0.75;
  testrel(22, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,1);
  f := 0.9375;
  testrel(23, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,2);
  f := 0.984375;
  testrel(24, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,4);
  f := 0.9990234375;
  testrel(25, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,6);
  f := 0.99993896484375;
  testrel(26, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,8);
  f := 0.999996185302734375;
  testrel(27, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,10);
  f := 0.9999997615814208984;
  testrel(28, NE, y, f, cnt,failed);

  p := 0.25;
  y := negbinom_cdf(p,r,1);
  f := 0.4375;
  testrel(29, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,2);
  f := 0.578125;
  testrel(30, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,5);
  f := 0.822021484375;
  testrel(31, NE, y, f, cnt,failed);

  y := negbinom_cdf(p,r,20);
  f := 0.9976215910457995051;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_hypergeo_pmf;
var
  n,n1,n2: longint;
  f,y: double;
  cnt, failed: integer;
const
  NE  = 1;
begin

  cnt := 0;
  failed := 0;
  writeln('Function: ','hypergeo_pmf');

  {Maple: stats[statevalf,pf,hypergeometric[n1, n2, n]](k);}

  n1 := 100;
  n2 := 200;
  n  := 50;

  y := hypergeo_pmf(n1,n2,n,0);
  f := 0.1458070693019822597e-9;
  testrel(1, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,5);
  f := 0.3329662375051472972e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,10);
  f := 0.1140119535468728544e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,15);
  f := 0.1145256680045404870;
  testrel(4, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,16);
  f := 0.1282808216541822096;
  testrel(5, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,17);
  f := 0.1290489702868419833;
  testrel(6, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,20);
  f := 0.7054335486599437926e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,30);
  f := 0.1522612415136588525e-4;
  testrel(8, NE, y, f, cnt,failed);


  n1 := 200;
  n2 := 100;
  n  := 50;

  y := hypergeo_pmf(n1,n2,n,0);
  f := 0.3241247045705749420e-28;
  testrel(9, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,10);
  f := 0.9914657919600567030e-13;
  testrel(10, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,20);
  f := 0.1522612415136588525e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,30);
  f := 0.7054335486599437926e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,33);
  f := 0.1290489702868419833;
  testrel(13, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,35);
  f := 0.1145256680045404870;
  testrel(14, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,40);
  f := 0.1140119535468728544e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,50);
  f := 0.1458070693019822597e-9;
  testrel(16, NE, y, f, cnt,failed);

  n1 := 200;
  n2 := 150;
  n  := 40;

  y := hypergeo_pmf(n1,n2,n,0);
  f := 0.6305576156685422360e-16;
  testrel(17, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,5);
  f := 0.6854951277117298538e-9;
  testrel(18, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,10);
  f := 0.1033878108095313435e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,20);
  f := 0.8380353022803073629e-1;
  testrel(20, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,21);
  f := 0.1096665608218504728;
  testrel(21, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,22);
  f := 0.1284352525327525682;
  testrel(22, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,23);
  f := 0.1345232262552923270;
  testrel(23, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,24);
  f := 0.1258645484459498172;
  testrel(24, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,25);
  f := 0.1050176499033465734;
  testrel(25, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,30);
  f := 0.6852689231387408111e-2;
  testrel(26, NE, y, f, cnt,failed);

  y := hypergeo_pmf(n1,n2,n,40);
  f := 0.2932117760869684499e-10;
  testrel(27, NE, y, f, cnt,failed);

  {Larger values from Boost test data}
  y := hypergeo_pmf(2598,600,2284,1692);
  f := 0.1612952941179119320e-83;
  testrel(28, NE, y, f, cnt,failed);

  y := hypergeo_pmf(2598,600,2284,1812);
  f := 0.2076831972632489092e-5;
  testrel(29, NE, y, f, cnt,failed);

  y := hypergeo_pmf(2598,600,2284,1848);
  f := 0.3043400366814929190e-1;
  testrel(30, NE, y, f, cnt,failed);

  y := hypergeo_pmf(2598,600,2284,1870);
  f := 0.1379967106042925267e-1;
  testrel(31, NE, y, f, cnt,failed);

  y := hypergeo_pmf(30197,719,24162,23557);
  f := 0.9957825121836800402e-5;
  testrel(32, NE, y, f, cnt,failed);

  y := hypergeo_pmf(30197,719,24162,23574);
  f := 0.2023261955389159984e-2;
  testrel(33, NE, y, f, cnt,failed);

  y := hypergeo_pmf(244,1029923,523360,179);
  f := 0.3252851387472271981e-12;
  testrel(34, NE, y, f, cnt,failed);

  y := hypergeo_pmf(244,2063033,523360,7);
  f := 0.4856769030435700293e-21;
  testrel(35, NE, y, f, cnt,failed);

  y := hypergeo_pmf(312,1029855,523360,149);
  f := 0.2530942714690769707e-1;
  testrel(36, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hypergeo_cdf;
var
  n,n1,n2: longint;
  f,y: double;
  cnt, failed: integer;
const
  NE  = 1;
begin

  cnt := 0;
  failed := 0;
  writeln('Function: ','hypergeo_cdf');

  {Maple: sum(stats[statevalf,pf,hypergeometric[n1, n2, n]](j), j=0..k);}
  n1 := 100;
  n2 := 200;
  n  := 50;

  y := hypergeo_cdf(n1,n2,n,0);
  f := 0.1458070693019822597e-9;
  testrel(1, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,5);
  f := 0.4001169492144623170e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,10);
  f := 0.1886559586767971603e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,15);
  f := 0.3548034277806054734;
  testrel(4, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,16);
  f := 0.4830842494347876830;
  testrel(5, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,17);
  f := 0.6121332197216296663;
  testrel(6, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,20);
  f := 0.8950818719479068637;
  testrel(7, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,30);
  f := 0.9999951348562412553;
  testrel(8, NE, y, f, cnt,failed);

  n1 := 200;
  n2 := 100;
  n  := 50;

  y := hypergeo_cdf(n1,n2,n,0);
  f := 0.3241247045705749420e-28;
  testrel(9, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,10);
  f := 0.1072725879761216243e-12;
  testrel(10, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,20);
  f := 0.2009126791011058333e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,30);
  f := 0.1754614829180875156;
  testrel(12, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,33);
  f := 0.5169157505652123170;
  testrel(13, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,35);
  f := 0.7597222402239350136;
  testrel(14, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,40);
  f := 0.9925355994870075694;
  testrel(15, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,45);
  f := 0.9999932849288290685;
  testrel(16, NE, y, f, cnt,failed);

  n1 := 200;
  n2 := 150;
  n  := 40;

  y := hypergeo_cdf(n1,n2,n,0);
  f := 0.6305576156685422360e-16;
  testrel(17, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,5);
  f := 0.7450135325309506333e-9;
  testrel(18, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,10);
  f := 0.128609830481474961e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,20);
  f := 0.2112666145268263496;
  testrel(20, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,21);
  f := 0.3209331753486768224;
  testrel(21, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,22);
  f := 0.4493684278814293906;
  testrel(22, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,23);
  f := 0.5838916541367217176;
  testrel(23, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,24);
  f := 0.7097562025826715348;
  testrel(24, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,25);
  f := 0.8147738524860181082;
  testrel(25, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,30);
  f := 0.9961139248041637270;
  testrel(26, NE, y, f, cnt,failed);

  y := hypergeo_cdf(n1,n2,n,35);
  f := 0.9999978850759868998;
  testrel(27, NE, y, f, cnt,failed);

  {Larger values from Boost test data}
  y := hypergeo_cdf(2598,600,2284,1692);
  f := 1.654453685667568298e-84;
  testrel(28, NE, y, f, cnt,failed);

  y := hypergeo_cdf(2598,600,2284,1812);
  f := 5.354493333876393888e-6;
  testrel(29, NE, y, f, cnt,failed);

  y := hypergeo_cdf(2598,600,2284,1848);
  f := 0.2426930943591920128;
  testrel(30, NE, y, f, cnt,failed);

  y := hypergeo_cdf(2598,600,2284,1870);
  f := 0.9332289351942023402;
  testrel(31, NE, y, f, cnt,failed);

  y := hypergeo_cdf(30197,719,24162,23557);
  f := 0.28832290318908713238e-4;
  testrel(33, NE, y, f, cnt,failed);

  y := hypergeo_cdf(30197,719,24162,23574);
  f := 0.0087632967578790970761218327;
  testrel(34, NE, y, f, cnt,failed);

  y := hypergeo_cdf(244,1029923,523360,179);
  f := 0.9999999999998102642;
  testrel(35, NE, y, f, cnt,failed);

  y := hypergeo_cdf(244,2063033,523360,7);
  f := 0.53100513491698108966e-21;
  testrel(36, NE, y, f, cnt,failed);

  y := hypergeo_cdf(312,1029855,523360,149);
  f := 0.1538437456126160276;
  testrel(37, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_rayleigh_pdf;
var
  y,f,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','rayleigh_pdf');

  b := 0.75;
  y := rayleigh_pdf(b, 1e-5);
  f := 0.1777777777619753086e-4;
  testrel(1, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 0.125);
  f := 0.2191571370542036040;
  testrel(2, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 0.75);
  f := 0.8087075462835112315;
  testrel(3, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 1.0);
  f := 0.7308662942349998861;
  testrel(4, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 2.0);
  f := 0.1015662250117346583;
  testrel(5, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 5.0);
  f := 0.1985450165513923870e-8;
  testrel(6, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 10.0);
  f := 0.4425104450360071126e-37;
  testrel(7, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 20.0);
  f := 0.1364871235187077375e-152;
  testrel(8, NE, y, f, cnt,failed);

  b := 2.5;
  y := rayleigh_pdf(b, 1e-5);
  f := 0.15999999999872e-5;
  testrel(9, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 0.125);
  f := 0.1997501561849161733e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 1.0);
  f := 0.1476986154218617253;
  testrel(11, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 2.0);
  f := 0.2323676918635810960;
  testrel(12, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 5.0);
  f := 0.1082682265892901535;
  testrel(13, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 10.0);
  f := 0.5367402046440189421e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 15.0);
  f := 0.3655195138731030825e-7;
  testrel(15, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 20.0);
  f := 0.4052532975710136231e-13;
  testrel(16, NE, y, f, cnt,failed);

  y := rayleigh_pdf(b, 30.0);
  f := 0.2582489356810146439e-30;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_rayleigh_cdf;
var
  y,f,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','rayleigh_cdf');

  b := 0.75;
  y := rayleigh_cdf(b, 1e-5);
  f := 0.8888888888493827161e-10;
  testrel(1, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 0.125);
  f := 0.1379288325608378218e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 0.75);
  f := 0.3934693402873665764;
  testrel(3, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 1.0);
  f := 0.5888877094928125641;
  testrel(4, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 2.0);
  f := 0.9714344992154496274;
  testrel(5, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 3.0);
  f := 0.9996645373720974882;
  testrel(6, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 5.0);
  f := 0.9999999997766368564;
  testrel(7, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 20.0);
  f := 1.0;
  testrel(8, NE, y, f, cnt,failed);

  b := 2.5;
  y := rayleigh_cdf(b, 1e-5);
  f := 0.7999999999968000000e-11;
  testrel(9, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 0.125);
  f := 0.1249219075419133499e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 1.0);
  f := 0.7688365361336421709e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 2.0);
  f := 0.2738509629263090751;
  testrel(12, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 5.0);
  f := 0.8646647167633873081;
  testrel(13, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 10.0);
  f := 0.9996645373720974882;
  testrel(14, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 15.0);
  f := 0.9999999847700202553;
  testrel(15, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 20.0);
  f := 0.9999999999999873358;
  testrel(16, NE, y, f, cnt,failed);

  y := rayleigh_cdf(b, 30.0);
  f := 1.0;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_rayleigh_inv;
var
  y,f,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','rayleigh_inv');

  b := 0.75;
  y := rayleigh_inv(b, 1e-8);
  f := 0.1060660174431471730e-3;
  testrel(1, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.125);
  f := 0.3875858830021908203;
  testrel(2, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.25);
  f := 0.5688957123306990994;
  testrel(3, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.5);
  f := 0.8830575168866060183;
  testrel(4, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.75);
  f := 1.248831916736546635;
  testrel(5, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.875);
  f := 1.529500485253213452;
  testrel(6, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.9);
  f := 1.609474519717010430;
  testrel(7, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.9990234375);
  f := 2.792473058294275600;
  testrel(8, NE, y, f, cnt,failed);

  b := 2.5;
  y := rayleigh_inv(b, 1e-8);
  f := 0.3535533914771572435e-3;
  testrel(9, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.125);
  f := 1.291952943340636068;
  testrel(10, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.25);
  f := 1.896319041102330331;
  testrel(11, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.5);
  f := 2.943525056288686728;
  testrel(12, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.75);
  f := 4.162773055788488782;
  testrel(13, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.875);
  f := 5.098334950844044839;
  testrel(14, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.9);
  f := 5.364915065723368099;
  testrel(15, NE, y, f, cnt,failed);

  y := rayleigh_inv(b, 0.9990234375);
  f := 9.308243527647585332;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{Maxwell Distribution

http://mathworld.wolfram.com/MaxwellDistribution.html
http://en.wikipedia.org/wiki/Maxwell-Boltzmann_distribution

Test values from Wolfram Alpha:

N[PDF[MaxwellDistribution[b],x],30]
pdf(1,0.5) = 0.176032663382149738887340220798
pdf(1,1.5) = 0.582829180496512774263448010796
pdf(1,3)   = 0.0797732714148841291608423485302
pdf(1,5)   = 0.0000743359757367148853954119816803

pdf(2,0.5) = 0.0241667573001780754338266095685
pdf(2,2)   = 0.241970724519143349797830192936
pdf(2,5)   = 0.109551878084803358513489513542
pdf(2,8)   = 0.00214128361223816562838519181505

pdf(5,0.5) = 0.00158781018990804706204211897400
pdf(5,6)   = 0.111851167670330653677364496083
pdf(5,10)  = 0.0863855464211008831209027206571
pdf(5,20)  = 0.000856513444895266251354076726020


N[CDF[MaxwellDistribution[b],x],30]
cdf(1,0.5)= 0.0308595957837267295007287796202
cdf(1,1.5)= 0.477832810464608685148713244176
cdf(1,2)  = 0.738535870050889377797177924024
cdf(1,5)  = 0.999984559501708898635097570099

cdf(2,0.5)= 0.00407859296442284501109470661429
cdf(2,1)  = 0.0308595957837267295007287796202
cdf(2,5)  = 0.899939166880605042855252180018
cdf(2,10) = 0.999984559501708898635097570099

cdf(5,0.5)= 0.000265165058655609828703287778194
cdf(5,5)  = 0.198748043098799197574804705393
cdf(5,10) = 0.738535870050889377797177924024
cdf(5,20) = 0.998866015710214677343299862579


N[InverseCDF[MaxwellDistribution[b],y],30]
Note: icdf(b,y) = b*icdf(1,y)!!

icdf(1,0.99)  = 3.36821417521872737743390060837
icdf(1,0.75)  = 2.02690526064547892102168810946
icdf(1,0.5)   = 1.53817225445505233444813196266
icdf(1,0.125) = 0.832080361133276623993541308689
icdf(1,0.01)  = 0.338868413841002651688751884974

icdf(2,0.99)  = 6.73642835043745475486780121674
icdf(2,0.75)  = 4.05381052129095784204337621892
icdf(2,0.5)   = 3.07634450891010466889626392531
icdf(2,0.125) = 1.66416072226655324798708261738
icdf(2,0.01)  = 0.677736827682005303377503769948

icdf(2,0.9921875) = 6.89293353359828505965776171799
icdf(2,0.875) = 4.79141457780476534844457091433
icdf(2,0.375) = 2.64875898849706732048133226452
icdf(2,0.250) = 2.20230143535862871466774425670
icdf(2,0.1/64)= 0.78965428521081280803013284852

icdf(5,0.99)  = 16.8410708760936368871695030418
icdf(5,0.75)  = 10.1345263032273946051084405473
icdf(5,0.5)   = 7.69086127227526167224065981328
icdf(5,0.125) = 4.16040180566638311996770654344
icdf(5,0.01)  = 1.69434206920501325844375942487
}


{---------------------------------------------------------------------------}
procedure test_maxwell_pdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Test maxwell_pdf');

  y := maxwell_pdf(1,0.5);
  f := 0.1760326633821497389;
  testrel(1, NE, y, f, cnt,failed);

  y := maxwell_pdf(1,1.5);
  f := 0.5828291804965127743;
  testrel(2, NE, y, f, cnt,failed);

  y := maxwell_pdf(1,3);
  f := 0.7977327141488412916e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := maxwell_pdf(1,5);
  f := 0.7433597573671488540e-4;
  testrel(4, NE, y, f, cnt,failed);

  y := maxwell_pdf(2,0.5);
  f := 0.2416675730017807543e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := maxwell_pdf(2,2);
  f := 0.2419707245191433498;
  testrel(6, NE, y, f, cnt,failed);

  y := maxwell_pdf(2,5);
  f := 0.1095518780848033585;
  testrel(7, NE, y, f, cnt,failed);

  y := maxwell_pdf(2,8);
  f := 0.2141283612238165628e-2;
  testrel(8, NE, y, f, cnt,failed);

  y := maxwell_pdf(5,0.5);
  f := 0.1587810189908047062e-2;
  testrel(9, NE, y, f, cnt,failed);

  y := maxwell_pdf(5,6);
  f := 0.1118511676703306537;
  testrel(10, NE, y, f, cnt,failed);

  y := maxwell_pdf(5,10);
  f := 0.8638554642110088312e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := maxwell_pdf(5,20);
  f := 0.8565134448952662514e-3;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_maxwell_cdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Test maxwell_cdf');

  y := maxwell_cdf(1,0.5);
  f := 0.3085959578372672950e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := maxwell_cdf(1,1.5);
  f := 0.4778328104646086851;
  testrel(2, NE, y, f, cnt,failed);

  y := maxwell_cdf(1,2);
  f := 0.7385358700508893778;
  testrel(3, NE, y, f, cnt,failed);

  y := maxwell_cdf(1,5);
       {0.999984559501708898635097570099}
  f := 2*0.4999922797508544493;
  testrel(4, NE, y, f, cnt,failed);

  y := maxwell_cdf(2,0.5);
  f := 0.4078592964422845011e-2;
  testrel(5, NE, y, f, cnt,failed);

  y := maxwell_cdf(2,1);
  f := 0.3085959578372672950e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := maxwell_cdf(2,5);
  f := 0.8999391668806050429;
  testrel(7, NE, y, f, cnt,failed);

  y := maxwell_cdf(2,10);
       {0.999984559501708898635097570099}
  f := 2*0.4999922797508544493;
  testrel(8, NE, y, f, cnt,failed);

  y := maxwell_cdf(5,0.5);
  f := 0.2651650586556098287e-3;
  testrel(9, NE, y, f, cnt,failed);

  y := maxwell_cdf(5,5);
       {0.198748043098799197574804705393}
  f := 0.3974960861975983951 / 2;
  testrel(10, NE, y, f, cnt,failed);

  y := maxwell_cdf(5,10);
  f := 0.7385358700508893778;
  testrel(11, NE, y, f, cnt,failed);

  y := maxwell_cdf(5,20);
  {    0.998866015710214677343299862579}
  f := 0.4994330078551073387*2;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_maxwell_inv;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 3;
  NE2 = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Test maxwell_inv');

  y := maxwell_inv(1,0.99);
  f := 3.368214175218727377;
  testrel(1, NE, y, f, cnt,failed);

  y := maxwell_inv(1,0.75);
  f := 2.026905260645478921;
  testrel(2, NE, y, f, cnt,failed);

  y := maxwell_inv(1,0.5);
  f := 1.538172254455052334;
  testrel(3, NE, y, f, cnt,failed);

  y := maxwell_inv(1,0.125);
  f := 0.8320803611332766240;
  testrel(4, NE, y, f, cnt,failed);

  y := maxwell_inv(1,0.01);
  f := 0.33886841384100265169;
  testrel(5, NE, y, f, cnt,failed);

  y := maxwell_inv(2,0.9921875);
  f := 6.89293353359828506;
  testrel(6, NE2, y, f, cnt,failed);

  y := maxwell_inv(2,0.875);
  f := 4.791414577804765348;
  testrel(7, NE, y, f, cnt,failed);

  y := maxwell_inv(2,0.375);
  f := 2.648758988497067320;
  testrel(8, NE, y, f, cnt,failed);

  y := maxwell_inv(2,0.25);
  f := 2.202301435358628715;
  testrel(9, NE, y, f, cnt,failed);

  y := maxwell_inv(2,1/64);
  f := 0.7896542852108128080;
  testrel(10, NE, y, f, cnt,failed);

  y := maxwell_inv(5,0.99);
  f := 16.84107087609363689;
  testrel(11, NE, y, f, cnt,failed);

  y := maxwell_inv(5,0.75);
       {10.1345263032273946051084405473;}
  f := 2*5.067263151613697303;
  testrel(12, NE, y, f, cnt,failed);

  y := maxwell_inv(5,0.5);
  f := 7.690861272275261672;
  testrel(13, NE, y, f, cnt,failed);

  y := maxwell_inv(5,0.125);
  f := 4.160401805666383120;
  testrel(14, NE, y, f, cnt,failed);

  y := maxwell_inv(5,0.01);
  f := 1.694342069205013258;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{
N[PDF[ExtremeValueDistribution[a,b],x],30]
N[CDF[ExtremeValueDistribution[a,b],x],30]
N[InverseCDF[ExtremeValueDistribution[a,b],y],30]

cdf(1,2,3) =     0.692200627555346353865421997183
pdf(1,2,3) =     0.127323190021791247909683347924
icdf(1,2,0.99) = 10.2002984535531599954458928018
}

{---------------------------------------------------------------------------}
procedure test_evt1_pdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Test evt1_pdf');
  {boost}
  y := evt1_pdf(1,3,5);
  f := 0.6750573240998512091e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := evt1_pdf(1,3,0);
  f := 0.1152223682858345643;
  testrel(2, NE, y, f, cnt,failed);

  y := evt1_pdf(0.5, 2, 0.125);
  f := 0.1805265483089020598;
  testrel(3, NE, y, f, cnt,failed);

  a := 1;
  b := 2;

  y := evt1_pdf(a,b,-1);
  f := 0.8968703936700859098e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,0);
  f := 0.158520960538971088;
  testrel(5, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,3);
  f := 0.1273231900217912479;
  testrel(6, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,10);
  f := 0.5493134841201401214e-2;
  testrel(7, NE, y, f, cnt,failed);

  a := 0;
  b := 1;

  y := evt1_pdf(a,b,0);
  f := 0.3678794411714423216;
  testrel(8, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,-2.5);
  f := 0.6236577187661991826e-4;
  testrel(9, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,2);
  f := 0.1182049515931431460;
  testrel(10, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,10);
  f := 0.4539786865564981977e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,700);
  f := 0.9859676543759770857e-304;
  testrel(12, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,-6.5);
  f := 0.9027618899817046825e-286;
  testrel(13, NE, y, f, cnt,failed);

  y := evt1_pdf(a,b,11000);
  f := 0;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_evt1_cdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Test evt1_cdf');

  {boost}
  y := evt1_cdf(0.5, 1.5, 0.125);
  f := 0.2769203340999089162;
  testrel(1, NE, y, f, cnt,failed);

  y := evt1_cdf(0.5, 2, -5);
  f := 1.608760113988777641e-7;
  testrel(2, NE, y, f, cnt,failed);

  y := evt1_cdf(0.5, 0.25, 0.75);
  f := 0.6922006275553463539;
  testrel(3, NE, y, f, cnt,failed);

  y := evt1_cdf(0.5, 0.25, 5);
       {0.999999984770020371326351248727041}
  f := 2*0.4999999923850101856;
  testrel(4, NE, y, f, cnt,failed);

  a := 1;
  b := 2;
  y := evt1_cdf(a,b,-5);
  f := 0.1892178694838292634e-8;
  testrel(5, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,0);
  f := 0.1922956455479649281;
  testrel(6, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,10);
  f := 0.9889524805037951500;
  testrel(7, NE, y, f, cnt,failed);

  a := 0;
  b := 1;

  y := evt1_cdf(a,b,-2);
  f := 0.6179789893310934986e-3;
  testrel(8, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,0);
  f := 0.3678794411714423216;
  testrel(9, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,2);
  f := 0.8734230184931166430;
  testrel(10, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,10);
        {0.999954601100798730506475075534}
  f := 2*0.4999773005503993653;
  testrel(11, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,20);
        {0.999999997938846379685619298220}
  f := 2*0.4999999989694231898;
  testrel(12, NE, y, f, cnt,failed);

  y := evt1_cdf(0,1,12000);
  f := 1.0;
  testrel(13, NE, y, f, cnt,failed);

  y := evt1_cdf(a,b,-8);
  f := 0.0;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_evt1_inv;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Test evt1_inv');

  a := 0;
  b := 1;
  y := evt1_inv(a,b,0.99);
  f := 4.600149226776579998;
  testrel(1, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.75);
  f := 1.245899323707238198;
  testrel(2, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.5);
  f := 0.3665129205816643270;
  testrel(3, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.125);
  f := -0.7320993680864453644;
  testrel(4, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.015625);
  f := -1.425246548646390674;
  testrel(5, NE, y, f, cnt,failed);

  {------------------------------------------------}
  a := 0.5;
  b := 0.25;
  y := evt1_inv(a,b,0.99);
  f := 1.650037306694145;
  testrel(6, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.75);
  f := 0.8114748309268095496;
  testrel(7, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.5);
  f := 0.5916282301454160818;
  testrel(8, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.125);
  f := 0.3169751579783886589;
  testrel(9, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.015625);
  f := 0.1436883628384023316;
  testrel(10, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,1e-300);
  f := -1.134453729976039215;
  testrel(11, NE, y, f, cnt,failed);

  {Alpha}
  a := 1;
  b := 2;

  y := evt1_inv(a,b,0.99);
  f := 10.20029845355316;
  testrel(12, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.75);
  f := 3.491798647414476397;
  testrel(13, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.5);
  f := 1.733025841163328654;
  testrel(14, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.125);
  f := -0.464198736172890728765612157380;
  testrel(15, NE, y, f, cnt,failed);

  y := evt1_inv(a,b,0.015625);
  f := -1.850493097292781348;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{Maple definitions:

kum_pdf := (a,b,x) -> a*b*x^(a-1)*(1-x^a)^(b-1);
kum_cdf := (a,b,x) -> 1-(1-x^a)^b;
kum_inv := (a,b,x) -> (1-(1-x)^(1/b))^(1/a);
}

{---------------------------------------------------------------------------}
procedure test_kumaraswamy_pdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Test kumaraswamy_pdf');

  a := 2;
  b := 5;

  y := kumaraswamy_pdf(a,b,0.0009765625);
  f := 0.9765587747150306052e-2;
  testrel(1, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.2);
  f := 1.698693120;
  testrel(2, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.4);
  f := 1.991485440;
  testrel(3, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.6);
  f := 1.006632960;
  testrel(4, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.8);
  f := 0.13436928;
  testrel(5, NE1, y, f, cnt,failed);  {???????}

  y := kumaraswamy_pdf(a,b,0.9990234375);
  f := 0.1450933120938096265e-9;
  testrel(6, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,1);
  f := 0.0;
  testrel(7, NE, y, f, cnt,failed);

  a := 0.5;
  b := 0.75;

  y := kumaraswamy_pdf(a,b,0.0009765625);
  f := 12.09562508943706705;
  testrel(8, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.2);
  f := 0.9724716100993977043;
  testrel(9, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.4);
  f := 0.7615068239595288480;
  testrel(10, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.6);
  f := 0.7026119956750337847;
  testrel(11, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.8);
  f := 0.7355263824425929538;
  testrel(12, NE, y, f, cnt,failed);

  y := kumaraswamy_pdf(a,b,0.9990234375);
  f := 2.523767830828681020;
  testrel(13, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_kumaraswamy_cdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Test kumaraswamy_cdf');

  a := 2;
  b := 5;

  y := kumaraswamy_cdf(a,b,0.0009765625);
  f := 0.4768362487092905884e-5;
  testrel(1, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.2);
  f := 0.1846273024;
  testrel(2, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.4);
  f := 0.5817880576;
  testrel(3, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.6);
  f := 0.8926258176;
  testrel(4, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.8);
  f := 0.9939533824;
  testrel(5, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.9990234375);
  f := 0.4999999999999858238*2;  {=0.9999999999999716476;}
  testrel(6, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,1);
  f := 1.0;
  testrel(7, NE, y, f, cnt,failed);

  a := 0.5;
  b := 0.75;

  y := kumaraswamy_cdf(a,b,0.0009765625);
  f := 0.2353026621732010763e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.2);
  f := 0.3589114579885191573;
  testrel(9, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.4);
  f := 0.5279560695712345950;
  testrel(10, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.6);
  f := 0.6728700954294068727;
  testrel(11, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.8);
  f := 0.8147901594656429391;
  testrel(12, NE, y, f, cnt,failed);

  y := kumaraswamy_cdf(a,b,0.9990234375);
  f := 0.9967146466464416880;
  testrel(13, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_kumaraswamy_inv;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Test kumaraswamy_inv');

  a := 2;
  b := 5;

  y := kumaraswamy_inv(a, b, 1.0/1024.0);
  f := 0.1397815576817549282e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.125);
  f := 0.1623355148599093934;
  testrel(2, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.5);
  f := 0.3597908235403952890;
  testrel(3, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.75);
  f := 0.4920789740933877362;
  testrel(4, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.875);
  f := 0.5833061328441120293;
  testrel(5, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.9);
  f := 0.6074888110243733127;
  testrel(6, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.99);
  f := 0.7758175232917227355;
  testrel(7, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.9990234375);
  f := 0.8660254037844386468;
  testrel(8, NE, y, f, cnt,failed);

  a := 0.5;
  b := 0.75;

  y := kumaraswamy_inv(a, b, 1.0/1024.0);
  f := 0.1694869037200626253e-5;
  testrel(9, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.125);
  f := 0.2659919863324392842e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.5);
  f := 0.3637896052527594082;
  testrel(11, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.75);
  f := 0.7098228789632848256;
  testrel(12, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.875);
  f := 0.87890625;
  testrel(13, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.9);
  f := 0.9093226580177763059;
  testrel(14, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.99);
  f := 0.4978478861043849227*2.0;  {=0.9956957722087698453}
  testrel(15, NE, y, f, cnt,failed);

  y := kumaraswamy_inv(a, b, 0.9990234375);
  f := 0.9998062348446667412;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{Maple definitions:
ks_cdf := x->exp(3/8*(Pi/x)^2)*exp(-(Pi/x)^2/2)^(3/4)*sqrt(Pi/2)/x*JacobiTheta2(0,exp(-(Pi/x)^2/2));
ks_inv := y->fsolve(ks_cdf(x)=y, x=0.0001..5);
c.f. http://stats.stackexchange.com/questions/45966/kolmogorov-smirnov-cdf
}

{---------------------------------------------------------------------------}
procedure test_kolmogorov_cdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
  NE1 = 4;
  NE2 = 64;
begin
  cnt := 0;
  failed := 0;
  writeln('Test kolmogorov_cdf');

  y := kolmogorov_cdf(0.01);
  f := 0.0;
  testabs(1, 0, y, f, cnt,failed);

  y := kolmogorov_cdf(0.1);
  f := 0.6609305242245470751e-52;
  testrel(2, NE2, y, f, cnt,failed);

  y := kolmogorov_cdf(0.125);
  f := 0.1027216723999461911e-32;
  testrel(3, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(0.375);
  f := 0.1035144856961865651e-2;
  testrel(4, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(0.4);
  f := 0.2807673222701733158e-2;
  testrel(5, NE1, y, f, cnt,failed);

  y := kolmogorov_cdf(0.425);
  f := 0.6373723693931488075e-2;
  testrel(6, NE1, y, f, cnt,failed);

  y := kolmogorov_cdf(0.45);
  f := 0.12589373847063302089e-1;
  testrel(7, NE1, y, f, cnt,failed);

  y := kolmogorov_cdf(0.475);
  f := 0.2226930807996159592e-1;
  testrel(8, NE1, y, f, cnt,failed);

  y := kolmogorov_cdf(0.5);
  f := 0.3605475633512490561e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1);
  f := 0.7300003283226454788;
  testrel(10, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.1796875);
  f := 0.8763641150488026113;
  testrel(11, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.25);
  f := 0.9121335860583089354;
  testrel(12, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.125);
  f := 0.8409611129029082319;
  testrel(13, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.5);
  f := 0.9777820373834748713;
  testrel(14, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.15);
  f := 0.8580401311219991918;
  testrel(15, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(1.2);
  f := 0.8877503333292750391;
  testrel(16, NE, y, f, cnt,failed);

  {16-bit compiler are inaccurate for constants 0.999....}
  y := kolmogorov_cdf(2);
 {$ifdef BIT16}
  f := 0.3331096915814067682*3.0;
 {$else}
  f := 0.9993290747442203047;
 {$endif}
  testrel(17, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(2.5);
 {$ifdef BIT16}
  f := 0.3333308488978852809*3.0;
 {$else}
  f := 0.9999925466936558427;
 {$endif}
  testrel(18, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(3);
 {$ifdef BIT16}
  f := 0.3333333231800135035*3;
 {$else}
  f := 0.9999999695400405106;
 {$endif}
  testrel(19, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(4);
 {$ifdef BIT16}
  f := 0.3333333333333248906*3.0;
 {$else}
  f := 0.9999999999999746717;
 {$endif}
  testrel(20, NE, y, f, cnt,failed);

  y := kolmogorov_cdf(4.75);
 {$ifdef BIT16}
  f := 0.3333333333333333333*3.0;
 {$else}
  f := 0.9999999999999999999;
 {$endif}
  testrel(21, NE, y, f, cnt,failed);

  {x > xmax, x <= xmin}
  y := kolmogorov_cdf(4.8);
  f := 1.0;
  testabs(22, 0, y, f, cnt,failed);

  y := kolmogorov_cdf(5);
  f := 1.0;
  testabs(23, 0, y, f, cnt,failed);

  y := kolmogorov_cdf(21/2048);
  f := 0.0;
  testabs(24, 0, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_kolmogorov_inv;
var
  x,y,f: double;
  n, cnt, failed: integer;
const
  NTEST = 500;
const
  NE = 1;
  NE1 = 4;  {tests cdf(inv(y))=y, y = 1/NTEST .. 1}
begin
  cnt := 0;
  failed := 0;
  writeln('Test kolmogorov_inv');

  x := 0.01;
  y := kolmogorov_inv(x);
  f := 0.4410276985179293684;
  testrel(1, NE, y, f, cnt,failed);

  x := 0.0625;
  y := kolmogorov_inv(x);
  f := 0.5345255069097582718;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.1;
  y := kolmogorov_inv(x);
  f := 0.5711732651063401632;
  testrel(3, NE, y, f, cnt,failed);

  x := 0.125;
  y := kolmogorov_inv(x);
  f := 0.5917613394822712151;
  testrel(4, NE, y, f, cnt,failed);

  x := 0.5;
  y := kolmogorov_inv(x);
  f := 0.8275735551899076901;
  testrel(5, NE, y, f, cnt,failed);

  x := 0.875;
  y := kolmogorov_inv(x);
  f := 1+0.1773581385807907705;
  testrel(6, NE, y, f, cnt,failed);

  x := 0.9;
  y := kolmogorov_inv(x);
  f := 1.223847870217082388;
  testrel(7, NE, y, f, cnt,failed);

  x := 0.99;
  y := kolmogorov_inv(x);
  f := 1.627623611518950347;
  testrel(8, NE, y, f, cnt,failed);

  for n:=1 to NTEST do begin
    {Avoid FPC single eval bug}
    x := NTEST;
    f := n/x;
    x := kolmogorov_inv(f);
    y := kolmogorov_cdf(x);
    testrel(10+n, NE1, y, f, cnt,failed);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.

{Part 3d of regression test for SPECFUN unit  (c) 2014-2015  W.Ehrhardt}

unit t_sfd3d;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_levy_cdf;
procedure test_levy_inv;
procedure test_levy_pdf;

procedure test_logseries_cdf;
procedure test_logseries_pmf;

procedure test_zipf_pmf;
procedure test_zipf_cdf;

procedure test_wald_cdf;
procedure test_wald_inv;
procedure test_wald_pdf;

procedure test_nakagami_cdf;
procedure test_nakagami_pdf;
procedure test_nakagami_inv;


implementation

uses
  specfun, t_sfd0;

{
N[PDF[LevyDistribution[1, 2], 3.5], 30] = 0.0956747327738255849972956318278
N[CDF[LevyDistribution[1, 2], 3.5], 30] = 0.371093369522697573790130103733
N[InverseCDF[LevyDistribution[1, 2], 0.875], 30] = 81.8190150737960697550460142716
}
{---------------------------------------------------------------------------}
procedure test_levy_pdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','levy_pdf');

  {levy_pdf := (m,s,x)-> sqrt(s/evalf(2*Pi))*exp(-s/(2*(x-m)))/(x-m)^(3/2);}
  a := -1/2;
  b := 2;

  y := levy_pdf(a, b, -0.25);
  f := 0.8266794141636821543e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 0);
  f := 0.2159638660527522078;
  testrel(2, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 0.5);
  f := 0.2075537487102973517;
  testrel(3, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 1);
{$ifdef BIT16}
  f := 0.3153468637585578479*0.5;
{$else}
  f := 0.15767343187927892396;
{$endif}
  testrel(4, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 2);
  f := 0.9567473277382558500e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 10);
  f := 0.1507577807792679211e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 100);
  f := 0.5544400996270593049e-3;
  testrel(7, NE, y, f, cnt,failed);

  a := 0;
  b := 1/2;

  y := levy_pdf(a, b, 0.125);
  f := 0.8638554642110088312;
  testrel(8, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 0.5);
  f := 0.4839414490382866996;
  testrel(9, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 1);
  f := 0.2196956447338611985;
  testrel(10, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 2);
  f := 0.8801633169107486944e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 5);
  f := 0.2400077896860271960e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 10);
  f := 0.8700369673862929858e-2;
  testrel(13, NE, y, f, cnt,failed);

  y := levy_pdf(a, b, 100);
  f := 0.2813904356065047971e-3;
  testrel(14, NE, y, f, cnt,failed);

  {Wolfram Alpha}
  y := levy_pdf(1,2,3.5);
  f := 0.9567473277382558500e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := levy_pdf(0,1,4);
  f := 0.4400816584553743472e-1;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_levy_cdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','levy_cdf');
  {levy_cdf := (m,s,x)->erfc(sqrt(s/(2*(x-m))));}

  a := -1/2;
  b := 2;
  y := levy_cdf(a, b, -0.25);
  f := 0.4677734981047265838e-2;
  testrel(1, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 0);
  f := 0.4550026389635841440e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 0.5);
  f := 0.1572992070502851307;
  testrel(3, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 1);
  f := 0.2482130789899235835;
  testrel(4, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 2);
  f := 0.3710933695226975738;
  testrel(5, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 10);
  f := 0.6625205835400574706;
  testrel(6, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 100);
  f := 0.8878153358231549984;
  testrel(7, NE, y, f, cnt,failed);

  a := 0;
  b := 1/2;

  y := levy_cdf(a, b, 0.125);
  f := 0.4550026389635841440e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 0.5);
  f := 0.3173105078629141028;
  testrel(9, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 1);
  f := 0.4795001221869534623;
  testrel(10, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 2);
  f := 0.6170750774519737927;
  testrel(11, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 5);
  f := 0.7518296340458492825;
  testrel(12, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 10);
  f := 0.8230632737581214761;
  testrel(13, NE, y, f, cnt,failed);

  y := levy_cdf(a, b, 100);
  f := 0.9436280222029833762;
  testrel(14, NE, y, f, cnt,failed);

  {Wolfram Alpha}
  y := levy_cdf(1,2,3.5);
  f := 0.3710933695226975738;
  testrel(15, NE, y, f, cnt,failed);

  y := levy_cdf(0,1,4);
  f := 0.6170750774519737927;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_levy_inv;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 3;
  NE1 = 9;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','levy_inv');
  {levy_inv := (m,s,x)-> m + s/((statevalf[icdf,normald[0,1]](x/2))^2);}

  a := -1/2;
  b := 2;

  y := levy_inv(a, b, 0.99);
  f := 12731.22877021246245;
  testrel(1, NE1, y, f, cnt,failed);

  y := levy_inv(a, b, 0.95);
  f := 508.1288891011180108;
  testrel(2, NE1, y, f, cnt,failed);

  y := levy_inv(a, b, 0.9);
{$ifdef BIT16}
  f := 252.3124707080669755*0.5;
{$else}
  f := 126.15623535403348776;
{$endif}
  testrel(3, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.75);
  f := 19.19840864364874886;
  testrel(4, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.5);
  f := 3.896218676635464808;
  testrel(5, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.125);
  f := 0.3497888411062026961;
  testrel(6, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 1e-3);
  f := -0.3152862823947457571;
  testrel(7, NE, y, f, cnt,failed);


  a := 0;
  b := 1/2;

  y := levy_inv(a, b, 0.99);
  f := 3182.932192553115612;
  testrel(8, NE1, y, f, cnt,failed);

  y := levy_inv(a, b, 0.95);
  f := 127.1572222752795027;
  testrel(9, NE1, y, f, cnt,failed);

  y := levy_inv(a, b, 0.9);
  f := 31.66405883850837194;
  testrel(10, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.75);
  f := 4.924602160912187214;
  testrel(11, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.5);
  f := 1.099054669158866202;
  testrel(12, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 0.125);
  f := 0.2124472102765506740;
  testrel(13, NE, y, f, cnt,failed);

  y := levy_inv(a, b, 1e-3);
  f := 0.4617842940131356072e-1;
  testrel(14, NE, y, f, cnt,failed);

  {Wolfram Alpha}
  y := levy_inv(1, 2, 0.875);
  f := 81.81901507379606976;
  testrel(15, NE, y, f, cnt,failed);

  y := levy_inv(0, 1, 0.875);
  f := 40.40950753689803488;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_logseries_cdf;
var
  a,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','logseries_cdf');

  a := 0.5;
  y := logseries_cdf(a,2);
  f := 0.9016844005556021296;
  testrel(1, NE, y, f, cnt,failed);

  y := logseries_cdf(a,5);
  f := 1-0.6644352054578320557e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := logseries_cdf(a,10);
  f := 0.9998812309831727881;
  testrel(3, NE, y, f, cnt,failed);

  y := logseries_cdf(a,20);
  f := 0.9999999372299021140;
  testrel(4, NE, y, f, cnt,failed);

  y := logseries_cdf(a,50);
  f := 1-0.24659029992140e-16;
  testrel(5, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.75;
  y := logseries_cdf(a,1);
  f := 0.5410106403333612778;
  testrel(6, NE, y, f, cnt,failed);

  y := logseries_cdf(a,2);
  f := 0.7438896304583717569;
  testrel(7, NE, y, f, cnt,failed);

  y := logseries_cdf(a,5);
  f := 0.9366246710771317121;
  testrel(8, NE, y, f, cnt,failed);

  y := logseries_cdf(a,10);
{$ifdef BIT16}
  f := 1 - 0.9115264471914481194e-2;
{$else}
  f := 0.9908847355280855188;
{$endif}
  testrel(9, NE, y, f, cnt,failed);

  y := logseries_cdf(a,100);
{$ifdef BIT16}
  f := 1.0 - 0.6680549048715623e-14;
{$else}
  f := 0.9999999999999933195;
{$endif}
  testrel(10, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.9990234375;
  y := logseries_cdf(a,1);
  f := 0.1441286159013095279;
  testrel(11, NE, y, f, cnt,failed);

  y := logseries_cdf(a,2);
  f := 0.2161225485512312306;
  testrel(12, NE, y, f, cnt,failed);

  y := logseries_cdf(a,5);
  f := 0.3287116142146085349;
  testrel(13, NE, y, f, cnt,failed);

  y := logseries_cdf(a,10);
  f := 0.4211550059363923490;
  testrel(14, NE, y, f, cnt,failed);

  y := logseries_cdf(a,100);
  f := 0.7346249735232929576;
  testrel(15, NE, y, f, cnt,failed);

  y := logseries_cdf(a,1000);
  f := 0.9671296153998283123;
  testrel(16, NE, y, f, cnt,failed);

  y := logseries_cdf(a,10000);
  f := 0.999999229589457260;
  testrel(17, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.015625;
  y := logseries_cdf(a,1);
  f := 0.992166994411624357;
  testrel(18, NE, y, f, cnt,failed);

  y := logseries_cdf(a,2);
  f := 0.999918299055465172;
  testrel(19, NE, y, f, cnt,failed);

  y := logseries_cdf(a,5);
  f := 0.999999999843904256;
  testrel(20, NE, y, f, cnt,failed);

  y := logseries_cdf(a,8);
  f := 0.999999999999999603;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_logseries_pmf;
var
  a,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
  NE1 = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','logseries_pmf');

  a := 0.5;
  y := logseries_pmf(a,2);
  f := 0.1803368801111204259;
  testrel(1, NE, y, f, cnt,failed);

  y := logseries_pmf(a,5);
  f := 0.9016844005556021296e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := logseries_pmf(a,10);
  f := 0.1408881875868128327e-3;
  testrel(3, NE, y, f, cnt,failed);

  y := logseries_pmf(a,20);
  f := 0.6879306034512345349e-7;
  testrel(4, NE, y, f, cnt,failed);

  y := logseries_pmf(a,50);
  f := 0.2562741203051934149e-16;
  testrel(5, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.75;
  y := logseries_pmf(a,1);
  f := 0.5410106403333612778;
  testrel(6, NE, y, f, cnt,failed);

  y := logseries_pmf(a,2);
  f := 0.2028789901250104792;
  testrel(7, NE, y, f, cnt,failed);

  y := logseries_pmf(a,5);
  f := 0.3423582958359551836e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := logseries_pmf(a,10);
  f := 0.4062161420319194805e-2;
  testrel(9, NE, y, f, cnt,failed);

  y := logseries_pmf(a,100);
  f := 0.2313507343989070722e-14;
  testrel(10, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.9990234375;
  y := logseries_pmf(a,1);
  f := 0.1441286159013095279;
  testrel(11, NE, y, f, cnt,failed);

  y := logseries_pmf(a,2);
  f := 0.7199393264992170266e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := logseries_pmf(a,5);
  f := 0.2871328753384213056e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := logseries_pmf(a,10);
  f := 0.1428667974926080843e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := logseries_pmf(a,100);
  f := 0.1308405153950717687e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := logseries_pmf(a,1000);
  f := 0.5430647467290243868e-4;       {????}
  testrel(16, NE1, y, f, cnt,failed);

  y := logseries_pmf(a,10000);
  f := 0.8240342489239981980e-9;
  testrel(17, NE, y, f, cnt,failed);

  {----------------------------}
  a := 0.015625;
  y := logseries_pmf(a,1);
  f := 0.992166994411624357;
  testrel(18, NE, y, f, cnt,failed);

  y := logseries_pmf(a,2);
  f := 0.7751304643840815288e-2;
  testrel(19, NE, y, f, cnt,failed);

  y := logseries_pmf(a,5);
  f := 0.1182755225195436903e-7;
  testrel(20, NE, y, f, cnt,failed);

  y := logseries_pmf(a,8);
  f := 0.2819908202160446413e-13;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{ Wolfram Alpha:
N[PDF[ZipfDistribution[2.5],7],30] = 0.000977992471032982553951402759293
N[CDF[ZipfDistribution[2.5],7],30] = 0.997710202865658627017645791782

N[PDF[ZipfDistribution[1.25],2],30] = 0.143968220424383466320360624049
N[CDF[ZipfDistribution[1.25],2],30] = 0.828800348676211005398401475338
}

{---------------------------------------------------------------------------}
procedure test_zipf_pmf;
var
  y,f,r: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zipf_pmf');

  y := zipf_pmf(2.5, 7);
  f := 0.0009779924710329825540; {Wolfram Alpha}
  testrel(1, NE, y, f, cnt,failed);

  y := zipf_pmf(1.25, 2);
  f := 0.1439682204243834663;    {Wolfram Alpha}
  testrel(2, NE, y, f, cnt,failed);

  r := 1.5;
  y := zipf_pmf(r, 1);
  f := 0.7454412962887771749;
  testrel(3, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 2);
  f := 0.1317766488955711757;
  testrel(4, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 5);
  {f:= 0.1333485929389813972e-1}
  f := 2*0.6667429646949069858e-2;
  testrel(5, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 10);
  f := 0.2357292358220957876e-2;
  testrel(6, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 100);
  f := 0.7454412962887771749e-5;
  testrel(7, NE, y, f, cnt,failed);

  r := 1.0;
  y := zipf_pmf(r, 1);
  f := 0.6079271018540266287;
  testrel(8, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 2);
  {f:= 0.1519817754635066572;}
  f := 2*0.7599088773175332858e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 5);
  f := 0.2431708407416106515e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 10);
  f := 0.6079271018540266287e-2;
  testrel(11, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 100);
  f := 0.6079271018540266287e-4;
  testrel(12, NE, y, f, cnt,failed);

  r := 0.5;
  y := zipf_pmf(r, 1);
  f := 0.3827933839994265622;
  testrel(13, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 2);
  f := 0.1353378988096702902;
  testrel(14, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 5);
  f := 0.3423808111839592446e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 10);
  f := 0.1210498966681642565e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := zipf_pmf(r, 100);
  f := 0.3827933839994265622e-3;
  testrel(17, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_zipf_cdf;
var
  y,f,r: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zipf_cdf');

  y := zipf_cdf(2.5, 7);
  f := 0.9977102028656586270;   {Wolfram Alpha}
  testrel(1, NE, y, f, cnt,failed);

  y := zipf_cdf(1.25, 2);
  f := 0.8288003486762110054;   {Wolfram Alpha}
  testrel(2, NE, y, f, cnt,failed);

  r := 1.5;
  y := zipf_cdf(r, 1);
  f := 0.7454412962887771749;
  testrel(3, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 2);
  f := 0.8772179451843483507;
  testrel(4, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 5);
  f := 0.961667926440314008;
  testrel(5, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 10);
  f := 0.985414381367689444;
  testrel(6, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 100);
  f := 0.999506750812669548;
  testrel(7, NE, y, f, cnt,failed);

  r := 1.0;
  y := zipf_cdf(r, 1);
  f := 0.6079271018540266287;
  testrel(8, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 2);
  f := 0.7599088773175332858;
  testrel(9, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 5);
  f := 0.8897688610191295296;
  testrel(10, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 10);
  f := 0.942145805354965341;
  testrel(11, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 100);
  f := 0.993951024017395072;
  testrel(12, NE, y, f, cnt,failed);

  r := 0.5;
  y := zipf_cdf(r, 1);
  f := 0.3827933839994265622;
  testrel(13, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 2);
  f := 0.5181312828090968525;
  testrel(14, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 5);
  f := 0.6738871580261133337;
  testrel(15, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 10);
  f := 0.7638016085053121605;
  testrel(16, NE, y, f, cnt,failed);

  y := zipf_cdf(r, 100);
  f := 0.923632241407361995;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{Maple:
wa_pdf := (mu,b,x) -> sqrt(b/(2*Pi*x^3))*exp(-b*(x-mu)^2/(2*mu^2*x));
wa_cdf := (mu,b,x) -> statevalf[cdf,normald](sqrt(b/x)*(x/mu-1))+exp(2.0*b/mu)*statevalf[cdf,normald](-sqrt(b/x)*(x/mu+1));
wa_inv := (mu,b,y) -> fsolve(evalf(wa_cdf(mu,b,x))=y, x=sqrt(1+2.25*(mu/b)^2)-mu/b);
}

{---------------------------------------------------------------------------}
procedure test_wald_cdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wald_cdf');
  {-------------------------------------}
  a := 1;
  b := 1;

  y := wald_cdf(a,b,0.0009765625);
  f := 0.2962614197217557007e-223;
  testrel(1, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.015625);
  f := 0.3356507274327110511e-14;
  testrel(2, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.125);
  f := 0.1206821184832047165e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.25);
  f := 0.1126907667166024006;
  testrel(4, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.5);
  f := 0.3649755481729598906;
  testrel(5, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,1);
  f := 0.6681020012231706064;
  testrel(6, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,2);
  f := 0.8854754259860064283;
  testrel(7, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,5);
  f := 0.990115297399673597;
  testrel(8, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,10);
  f := 0.999649585462791181;
  testrel(9, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,20);
  f := 0.999999055200388898;
  testrel(10, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 1;
  b := 0.25;

  y := wald_cdf(a,b,0.0009765625);
  f := 0.1640465876801375975e-56;
  testrel(11, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.015625);
  f := 0.8119013171481177982e-4;
  testrel(12, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.09375);
  f := 0.1305229581458689291;
  testrel(13, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.125);
  f := 0.1999708176969944981;
  testrel(14, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.25);
  f := 0.4008143814660667369;
  testrel(15, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.5);
  f := 0.5999487302745564651;
  testrel(16, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,1);
  f := 0.7615782918651233717;
  testrel(17, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,2);
  f := 0.8762751204427933981;
  testrel(18, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,5);
  f := 0.9626012216974486591;
  testrel(19, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,10);
  f := 0.990225131369742161;
  testrel(20, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,50);
  f := 0.9999908061013105845;
  testrel(21, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 2;
  b := 4;
  y := wald_cdf(a,b,0.00390625);
  f := 0.8041455158779551833e-223;
  testrel(22, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.015625);
  f := 0.9368460534686932296e-56;
  testrel(23, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.125);
  f := 0.107382590485147250511e-6;
  testrel(24, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.25);
  f := 0.4181357460631188740e-3;
  testrel(25, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,0.5);
  f := 0.2805684041471993641e-1;
  testrel(26, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,1);
  f := 0.2323571891918430400;
  testrel(27, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,2);
  f := 0.6276978381552528719;
  testrel(28, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,3);
  f := 0.8244079562051371025;
  testrel(29, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,5);
  f := 0.9577838788517624366;
  testrel(30, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,10);
  f := 0.9983288481644421095;
  testrel(31, NE, y, f, cnt,failed);

  y := wald_cdf(a,b,20);
  f := 0.9999952073513724241;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_wald_inv;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wald_inv');
  {-------------------------------------}
  a := 1;
  b := 1;

  y := wald_inv(a,b,1e-10);
  f := 0.2285389221845809094e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := wald_inv(a,b,1e-3);
  f := 0.7921847779047665021e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.125);
  f := 0.2617856914896014894;
  testrel(3, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.5);
  f := 0.6758413056952391192;
  testrel(4, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.75);
  f := 1.244059755876046385;
  testrel(5, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.875);
  f := 1.909444319675166507;
  testrel(6, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9);
  f := 2.143033912957148442;
  testrel(7, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9921875);      {Needs restart}
  f := 5.325414747748768454;
  testrel(8, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9990234375);   {Needs restart}
  f := 8.39143866628992312;
  testrel(9, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 1;
  b := 0.25;

  y := wald_inv(a,b,1e-10);
  f := 0.5908915560393917305e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := wald_inv(a,b,1e-3);
  f := 0.2215015383550524002e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.125);
  f := 0.9137917616427910634e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.5);
  f := 0.3491836091137611914;
  testrel(13, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.75);
  f := 0.9443802870883098166;
  testrel(14, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.875);
  f := 1.980908162655757539;
  testrel(15, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9);
  f := 2.422032013286069764;
  testrel(16, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9921875);      {Needs restart}
  f := 10.97342024988035286;
  testrel(17, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9990234375);   {Needs restart}
  f := 21.38342463473598503;
  testrel(18, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 2;
  b := 4;
  y := wald_inv(a,b,1e-10);
  f := 0.8761694487689471823e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := wald_inv(a,b,1e-3);
  f := 0.2791099375319904766;
  testrel(20, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.125);
  f := 0.7721265466305410259;
  testrel(21, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.5);
  f := 1.608678082592003238;
  testrel(22, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.75);
  f := 2.526706448319834535;
  testrel(23, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.875);
  f := 3.463215256979898856;
  testrel(24, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9);
  f := 3.771990672253054854;
  testrel(25, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9921875);      {Needs restart}
  f := 7.544464127793504471;
  testrel(26, NE, y, f, cnt,failed);

  y := wald_inv(a,b,0.9990234375);   {Needs restart}
  f := 10.87830366087049127;
  testrel(27, NE, y, f, cnt,failed);

  y := wald_inv(2,8,1e-10);          {Needs restart for samll y}
  f := 0.1621063575982417337;
  testrel(28, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_wald_pdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wald_pdf');
  {-------------------------------------}
  a := 1;
  b := 1;

  y := wald_pdf(a,b,0.0009765625);
  f := 0.1554775506885176576e-217;
  testrel(1, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.015625);
  f := 0.6976830341002265365e-11;
  testrel(2, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.125);
  f := 0.4221999674412029973;
  testrel(3, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.25);
  f := 1.036140765327133821;
  testrel(4, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.5);
  f := 0.8787825789354447941;
  testrel(5, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,1);
  f := 0.3989422804014326779;
  testrel(6, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,2);
  f := 0.1098478223669305993;
  testrel(7, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,5);
  f := 0.7204168934430732589e-2;
  testrel(8, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,10);
  f := 0.2197948003186267005e-3;
  testrel(9, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,100);
  f := 0.2081176820202829743e-24;
  testrel(10, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 1;
  b := 0.25;

  y := wald_pdf(a,b,0.0009765625);
  f := 0.2158524227137894615e-51;
  testrel(11, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.015625);
  f := 0.4390556446584921146e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.09375);
  f := 2.324599169820407270;
  testrel(13, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.125);
  f := 2.098980181160592875;
  testrel(14, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.25);
  f := 1.204549728619217620;
  testrel(15, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.5);
  f := 0.5300070646880571217;
  testrel(16, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,1);
  f := 0.1994711402007163390;
  testrel(17, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,2);
  f := 0.6625088308600714022e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,5);
  f := 0.1195934159672819812e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,10);
  f := 0.2291695475027150821e-2;
  testrel(20, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,1000);
  f := 0.4183975819237239642e-59;
  testrel(21, NE, y, f, cnt,failed);

  {-------------------------------------}
  a := 2;
  b := 4;
  y := wald_pdf(a,b,0.00390625);
  f := 0.1055032911078540971e-217;
  testrel(22, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.015625);
  f := 0.7703930869108508351e-52;
  testrel(23, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.125);
  f := 0.1410291505970886346e-4;
  testrel(24, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.25);
  f := 0.1396292312073216105e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,0.5);
  f := 0.2378605784472587431;
  testrel(26, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,1);
  f := 0.4839414490382866996;
  testrel(27, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,2);
  f := 0.2820947917738781435;
  testrel(28, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,5);
  f := 0.2901482939356917119e-1;
  testrel(29, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,10);
  f := 0.1028484425270353499e-2;
  testrel(30, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,50);
  f := 0.2225052137962172005e-12;
  testrel(31, NE, y, f, cnt,failed);

  y := wald_pdf(a,b,100);
  f := 0.1114600004544149534e-23;
  testrel(32, NE, y, f, cnt,failed);

  y := wald_pdf(2,8,0.125);
  f := 0.1557965141138566289e-10;
  testrel(33, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_nakagami_cdf;
var
  m,w,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','nakagami_cdf');

  m := 1/2;
  w := 1;

  y := nakagami_cdf(m,w,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,0.1);
  f := 0.7965567455405796293e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,0.5);
  f := 0.38292492254802620728;
  testrel(3, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1);
  f := 0.68268949213708589717;
  testrel(4, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1.5);
  f := 0.86638559746228386799;
  testrel(5, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2);
  f := 1.908999472207283171*0.5;
  testrel(6, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2.5);
  f := 1.975161338696895459*0.5;
  testrel(7, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,3);
  f := 1.994600407873479622*0.5;
  testrel(8, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,10);
  f := 1.0;
  testrel(9, NE, y, f, cnt,failed);


  {------------------------------------}
  m := 1;
  w := 1;

  y := nakagami_cdf(m,w,0);
  f := 0;
  testrel(10, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,0.1);
  f := 0.01990033250166389285*0.5;
  testrel(11, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,0.5);
  f := 0.2211992169285951318;
  testrel(12, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1);
  f := 0.6321205588285576784;
  testrel(13, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1.5);
  f := 0.8946007754381356632;
  testrel(14, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2);
  f := 1.963368722222531639*0.5;
  testrel(15, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2.5);
  f := 1.996139091727544582*0.5;
  testrel(16, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,3);
  f := 1.999753180391826641*0.5;
  testrel(17, NE, y, f, cnt,failed);

  {------------------------------------}
  m := 5;
  w := 2;

  y := nakagami_cdf(m,w,0);
  f := 0;
  testrel(18, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,0.1);
  f := 0.7970282082925807112e-10;
  testrel(19, NE+2, y, f, cnt,failed);  {???}

  y := nakagami_cdf(m,w,0.5);
  f := 0.0004739871032458523993;
  testrel(20, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1);
  f := 0.10882198108584875765;
  testrel(21, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,1.5);
  f := 0.6616246257772604542;
  testrel(22, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2);
  f := 1.941494623846077855*0.5;
  testrel(23, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,2.5);
  f := 1.998933083406074633*0.5;
  testrel(24, NE, y, f, cnt,failed);

  y := nakagami_cdf(m,w,3);
  f := 1.999995650545484666*0.5;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_nakagami_inv;
var
  m,w,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','nakagami_inv');

  m := 1/2;
  w := 1;

  y := nakagami_inv(m,w,1/8);
  f := 0.1573106846101706955;
  testrel(1, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/4);
  f := 0.3186393639643751630;
  testrel(2, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/2);
  f := 0.6744897501960817432;
  testrel(3, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,3/4);
  f := 1.150349380376008178;
  testrel(4, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,9/10);
  f := 1.644853626951472715;
  testrel(5, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,99/100);
  f := 2.575829303548900761;
  testrel(6, NE, y, f, cnt,failed);

  {------------------------------------}
  m := 1;
  w := 1;

  y := nakagami_inv(m,w,1/8);
  f := 0.3654194748840332718;
  testrel(7, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/4);
  f := 0.5363600213026516459;
  testrel(8, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/2);
  f := 0.8325546111576977564;
  testrel(9, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,3/4);
  f := 1.177410022515474691;
  testrel(10, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,9/10);
  f := 1.517427129385146351;
  testrel(11, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,99/100);
  f := 2.145966026289347240;
  testrel(12, NE, y, f, cnt,failed);

  {------------------------------------}
  m := 5;
  w := 2;

  y := nakagami_inv(m,w,1/100);
  f := 0.7152918509513730029;
  testrel(13, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/8);
  f := 1.0231438807461379470;
  testrel(14, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/4);
  f := 1.160792899009521174;
  testrel(15, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,1/2);
  f := 1.366880957917840586;
  testrel(16, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,3/4);
  f := 1.584226082154272387;
  testrel(17, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,9/10);
  f := 1.788137532300312792;
  testrel(18, NE, y, f, cnt,failed);

  y := nakagami_inv(m,w,99/100);
  f := 2.154495354321023347;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_nakagami_pdf;
var
  m,w,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','nakagami_pdf');

  m := 1/2;
  w := 1;

  y := nakagami_pdf(m,w,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1e-20);
  f := 0.7978845608028653559;
  testrel(2, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,0.5);
  f := 0.7041306535285989555;
  testrel(3, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1);
  f := 0.4839414490382866996;
  testrel(4, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1.5);
  f := 0.2590351913317834552;
  testrel(5, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2);
  f := 0.1079819330263761039;
  testrel(6, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2.5);
  f := 0.03505660098713707472;
  testrel(7, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,3);
  f := 0.008863696823876014351;
  testrel(8, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,10);
  f := 0.1538919725341283869e-21;
  testrel(9, NE, y, f, cnt,failed);

  {------------------------------------}
  m := 1;
  w := 1;

  y := nakagami_pdf(m,w,0);
  f := 0;
  testrel(10, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1e-20);
  f := 2e-20;
  testrel(11, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,0.5);
  f := 0.7788007830714048682;
  testrel(12, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1);
  f := 0.7357588823428846432;
  testrel(13, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1.5);
  f := 0.3161976736855930103;
  testrel(14, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2);
  f := 0.073262555554936721175;
  testrel(15, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2.5);
  f := 0.009652270681138546211;
  testrel(16, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,3);
  f := 0.0007404588245200772970;
  testrel(17, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,10);
  f := 0.7440151952041671926e-42;
  testrel(18, NE, y, f, cnt,failed);

  {------------------------------------}
  m := 5;
  w := 2;

  y := nakagami_pdf(m,w,0);
  f := 0;
  testrel(19, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1e-20);
  f := 0.8138020833333333333e-179;
  testrel(20, NE+1, y, f, cnt,failed);  {!!!}

  y := nakagami_pdf(m,w,0.5);
  f := 0.008507751282358014456;
  testrel(21, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1);
  f := 0.6680094289054263930;
  testrel(22, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,1.5);
  f := 1.128323590059242638;
  testrel(23, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2);
  f := 0.1891663740103535481;
  testrel(24, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,2.5);
  f := 0.005083087616261134155;
  testrel(25, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,3);
  f := 0.00002710093327704558179;
  testrel(26, NE, y, f, cnt,failed);

  y := nakagami_pdf(m,w,10);
  f := 0.2172192558220439773e-98;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.

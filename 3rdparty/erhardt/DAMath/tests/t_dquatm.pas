{Test/dev program for DAMQuat  (c) W.Ehrhardt 2017}
unit t_dquatm;

interface

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

{$x+}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


procedure test_quat;

implementation

uses
  damath, damquat;


type
  tq = array[0..3] of double;

var
  eps, eps_f, maxerr, maxerr_f: double;

var
  cnt, failed: integer;


{---------------------------------------------------------------------------}
function reldevq(const a: quaternion; const b: TQ): double;
var
  t: quaternion;
  r: double;
begin
  inc(cnt);
  t.r := b[0];
  t.x := b[1];
  t.y := b[2];
  t.z := b[3];
  r := qabs(t);
  qsub(a,t,t);
  if r<0.5 then r := qabs(t)
  else r := qabs(t)/r;
  reldevq := r;
  if r > eps then begin
    writeln(cnt:4, r);
    inc(failed);
  end;
  if r>maxerr then maxerr := r;
end;


{---------------------------------------------------------------------------}
function reldevx(x,y: double): double;
var
  r: double;
begin
  inc(cnt);
  r := abs(x-y);
  if r >= 0.5 then r := r/abs(y);
  reldevx := r;
  if r > eps then begin
    writeln(cnt:4, r);
    inc(failed);
  end;
  if r>maxerr then maxerr := r;
end;


var
  a,b,c: quaternion;

{Test case computed with Mathematica package Quaternions}

{---------------------------------------------------------------------------}
procedure test1;
var
  x,y: double;
const
  raddab: tq = (0, 2.25, 3.0, 3.5);
  rsubab: tq = (2.0, 1.75, 3.0, 4.5);
  rmulab: tq = (0.5, -3.25, -1.0, -5.25);
  rmulba: tq = (0.5, -0.25, -5.0, -3.75);
  rdivab: tq = (-1.9047619047619047619, -0.57142857142857142857, -3.8095238095238095238, -2.0952380952380952381);
  rconj:  tq = (1,-2,-3,-4);
  rneg:   tq = (-1,-2,-3,-4);
  racos:  tq = (1.3901159465978957317, -0.8918751005821986307, -1.3378126508732979461, -1.7837502011643972614);
  racosh: tq = (2.4014472020074010366, 0.5162761016176176326, 0.7744141524264264490, 1.0325522032352352653);
  racot:  tq = (0.034428244650522299737, -0.06731789442942041199, -0.10097684164413061798, -0.13463578885884082397);
  racoth: tq = (0.032302932870001552068, -0.06603374190320327607, -0.09905061285480491411, -0.13206748380640655214);
  racsc:  tq = (0.032814272582458290012, -0.06634907934692852819, -0.09952361902039279228, -0.13269815869385705638);
  racsch: tq = (0.033876580580246358868, -0.06699111776677093971, -0.10048667665015640956, -0.13398223553354187941);
  rasec:  tq = (1.5379820542124383292, 0.06634907934692852819, 0.09952361902039279228, 0.13269815869385705638);
  rasech: tq = (0.17865036354242713609, -0.5711921953344312951, -0.8567882930016469426, -1.1423843906688625901);
  rasin:  tq = (0.18068038019700088753, 0.8918751005821986307, 1.3378126508732979461, 1.7837502011643972614);
  rasinh: tq = (2.3858899025852419682, 0.5140526006627880378, 0.7710789009941820566, 1.0281052013255760755);
  ratan:  tq = (1.5363680821443743195, 0.06731789442942041199, 0.10097684164413061798, 0.13463578885884082397);
  ratanh: tq = (0.032302932870001552068, 0.5173453683196951253, 0.7760180524795426879, 1.0346907366393902505);
  rcbrt:  tq = (1.5776218631838213982, 0.2920419815730954907, 0.4380629723596432360, 0.5840839631461909813);
  rcos:   tq = (58.933646167943954133, -34.08618369046560641, -51.12927553569840962, -68.17236738093121282);
  rcosh:  tq = (0.96158511763695705371, -0.3413521745610165670, -0.5120282618415248505, -0.6827043491220331339);
  rcot:   tq = (0.000038214980491101106557, -0.3713841806354823630, -0.5570762709532235446, -0.7427683612709647261);
  rcoth:  tq = (0.91000464302431212833, 0.09083062976716011420, 0.13624594465074017130, 0.18166125953432022839);
  rcsc:   tq = (0.0077147758359599268279, -0.0018396437588678951504, -0.0027594656383018427257, -0.003679287517735790301);
  rcsch:  tq = (0.36749725482767571642, 0.22491676803326495510, 0.3373751520498974327, 0.4498335360665299102);
  rexp:   tq = (1.6939227236833002502, -0.78955962454155853119, -1.1843394368123377968, -1.5791192490831170624);
  rinv:   tq = (0.033333333333333333333, -0.066666666666666666667, -0.1, -0.13333333333333333333);
  rln:    tq = (1.7005986908310776877, 0.5151902926640850204, 0.7727854389961275306, 1.0303805853281700408);
  rpowx:  tq = (-66.50377063575604140, -8.36042820857835966, -12.54064231286753949, -16.72085641715671932);
  rsec:   tq = (0.0049537738746715510481, 0.0028651756209329915419, 0.004297763431399487313, 0.005730351241865983084);
  rsech:  tq = (0.54344484351170861139, 0.19291696146729491240, 0.28937544220094236860, 0.3858339229345898248);
  rsin:   tq = (91.783715784034691694,21.886486853029179758,32.829730279543769637,43.772973706058359516);
  rsinh:  tq = (0.73233760604634319647, -0.4482074499805419642, -0.6723111749708129463, -0.8964148999610839284);
  rsqrt:  tq = (1.7996146219471074767, 0.5556745248702424845, 0.8335117873053637267, 1.1113490497404849690);
  rtan:   tq = (0.000038216317250070778857, 0.3713971716439371360, 0.5570957574659057040, 0.7427943432878742720);
  rtanh:  tq = (1.024869536055662052, -0.10229568178876420056, -0.15344352268314630084, -0.20459136357752840112);
  runit:  tq = (0.18257418583505537115, 0.36514837167011074230, 0.54772255750516611346, 0.73029674334022148461);
begin
  with a do begin
    r := 1;
    x := 2;
    y := 3;
    z := 4;
  end;
  with b do begin
    r := -1;
    x := 0.25;
    y := 0;
    z := -0.5;
  end;

  x := qabs(a);
  y := 5.4772255750516611346;
  reldevx(x,y);

  x := qnorm(a);
  y := 30;
  reldevx(x,y);

  x := qabs_im(a);
  y := 5.3851648071345040313;
  reldevx(x,y);

  x := qarg(a);
  y := 1.3871923165159780480;
  reldevx(x,y);

  x := qdot(a,b);
  y := -2.5;
  reldevx(x,y);

  qadd(a,b,c);     reldevq(c, raddab);
  qsub(a,b,c);     reldevq(c, rsubab);
  qmul(a,b,c);     reldevq(c, rmulab);
  qmul(b,a,c);     reldevq(c, rmulba);
  qdiv(a,b,c);     reldevq(c, rdivab);

  qconj(a,c);      reldevq(c, rconj);
  qneg(a,c);       reldevq(c, rneg);
  qarccos(a,c);    reldevq(c, racos);
  qarccosh(a,c);   reldevq(c, racosh);
  qarccot(a,c);    reldevq(c, racot);
  qarccoth(a,c);   reldevq(c, racoth);
  qarccsc(a,c);    reldevq(c, racsc);
  qarccsch(a,c);   reldevq(c, racsch);
  qarcsec(a,c);    reldevq(c, rasec);
  qarcsech(a,c);   reldevq(c, rasech);
  qarcsin(a,c);    reldevq(c, rasin);
  qarcsinh(a,c);   reldevq(c, rasinh);
  qarctan(a,c);    reldevq(c, ratan);
  qarctanh(a,c);   reldevq(c, ratanh);
  qcbrt(a,c);      reldevq(c, rcbrt);
  qcos(a,c);       reldevq(c, rcos);
  qcosh(a,c);      reldevq(c, rcosh);
  qcot(a,c);       reldevq(c, rcot);
  qcoth(a,c);      reldevq(c, rcoth);
  qcsc(a,c);       reldevq(c, rcsc);
  qcsch(a,c);      reldevq(c, rcsch);
  qexp(a,c);       reldevq(c, rexp);
  qinv(a,c);       reldevq(c, rinv);
  qln(a,c);        reldevq(c, rln);
  qpowx(a,5/2,c);  reldevq(c, rpowx);
  qsec(a,c);       reldevq(c, rsec);
  qsech(a,c);      reldevq(c, rsech);
  qsin(a,c);       reldevq(c, rsin);
  qsinh(a,c);      reldevq(c, rsinh);
  qsqrt(a,c);      reldevq(c, rsqrt);
  qtan(a,c);       reldevq(c, rtan);
  qtanh(a,c);      reldevq(c, rtanh);
  qunit(a,c);      reldevq(c, runit);
end;


{---------------------------------------------------------------------------}
procedure test2;
const
  raddab: tq = (0.41421356237309504880, 0.25, 1.0, 2.6415926535897932385);
  rsubab: tq = (-2.4142135623730950488, 0.25, -1.0, -3.6415926535897932385);
  rmulab: tq = (0.15658276442180157043, 0.85355339059327376220, -1.7853981633974483096, -3.5986994347763407629);
  rmulba: tq = (0.15658276442180157043, -0.14644660940672623780, -0.21460183660255169038, -4.0986994347763407629);
  rdivab: tq = (-0.23194262979173804015, -0.011379262706344737364, 0.13872984030855849480, 0.16973993949793344448);
  rconj:  tq = (-1.0, -0.25, 0, 0.5);
  rneg:   tq = (1.0, -0.25, 0, 0.5);
  racos:  tq = (2.4322949749142515929, -0.34773308749246167828, 0, 0.6954661749849233566);
  racosh: tq = (0.77755482165902620194, 1.0877553810478824550, 0, -2.1755107620957649100);
  racot:  tq = (-0.70789979243547781846, -0.11795445269589311849, 0, 0.23590890539178623698);
  racoth: tq = (-0.65616714804078975190, -0.29029606105967261780, 0, 0.5805921221193452356);
  racsc:  tq = (-0.72395373211246824037, -0.24222374519305803231, 0, 0.4844474903861160646);
  racsch: tq = (-0.73853593745270040344, -0.15104661362092975032, 0, 0.30209322724185950064);
  rasec:  tq = (2.2947500589073648596, 0.24222374519305803231, 0, -0.4844474903861160646);
  rasech: tq = (0.54162876001626568082, -1.0262434246177029215, 0, 2.0524868492354058429);
  rasin:  tq = (-0.86149864811935497367, 0.34773308749246167828, 0, -0.6954661749849233566);
  rasinh: tq = (-0.9374497102973956420, 0.17414182849018498859, 0, -0.3482836569803699772);
  ratan:  tq = (-0.86289653435941880077, 0.11795445269589311849, 0, -0.23590890539178623698);
  ratanh: tq = (-0.65616714804078975190, 0.41218541204440002152, 0, -0.8243708240888000430);
  rcbrt:  tq = (0.6688803396845122365, 0.3598552142613759157, 0, -0.7197104285227518315);
  rcos:   tq = (0.62694606610827940714, 0.22149687669695598632, 0, -0.4429937533939119726);
  rcosh:  tq = (1.3081880562340542043, -0.27873552354684978129, 0, 0.5574710470936995626);
  rcot:   tq = (-0.43114585830418261574, -0.28964769313738117487, 0, 0.5792953862747623497);
  rcoth:  tq = (-1.090868916177308256, -0.12095766035088475325, 0, 0.24191532070176950649);
  rcsc:   tq = (-0.925935549219680973, -0.13486943380539163535, 0, 0.26973886761078327070);
  rcsch:  tq = (-0.59932929246568456603, -0.22016102584851973231, 0, 0.4403220516970394646);
  rexp:   tq = (0.31187967771506699914, 0.087254054192293704479, 0, -0.17450810838458740896);
  rinv:   tq = (-0.76190476190476190476, -0.19047619047619047619, 0, 0.38095238095238095238);
  rln:    tq = (0.13596685774182087942, 1.1770004316689132715, 0, -2.3540008633378265430);
  rpowx:  tq = (1.3435530163531774110, 0.1835298110671513810, 0, -0.367059622134302762);
  rsec:   tq = (0.98211113848730822670, -0.34697490183575079391, 0, 0.6939498036715015878);
  rsech:  tq = (0.62299906959716031065, 0.13274236147153439204, 0, -0.26548472294306878409);
  rsin:   tq = (-0.97641064629903734614, 0.14222150897963236432, 0, -0.28444301795926472863);
  rsinh:  tq = (-0.99630837851898720515, 0.36598957773914348577, 0, -0.7319791554782869715);
  rsqrt:  tq = (0.2698554462475790034, 0.4632109588231867385, 0, -0.9264219176463734770);
  rtan:   tq = (-0.71220730088217454856, 0.47846731625236014388, 0, -0.9569346325047202878);
  rtanh:  tq = (-0.86361079696450255571, 0.09575883949525825593, 0, -0.19151767899051651187);
  runit:  tq = (-0.87287156094396952506, 0.21821789023599238127, 0, -0.43643578047198476253);
var
  x,y: double;
begin
  with b do begin
    r := sqrt2;
    x := 0;
    y := 1;
    z := Pi;
  end;
  with a do begin
    r := -1;
    x := 0.25;
    y := 0;
    z := -0.5;
  end;

  x := qabs(a);
  y := 1.1456439237389600016;
  reldevx(x,y);

  x := qnorm(a);
  y := 1.3125;
  reldevx(x,y);

  x := qabs_im(a);
  y := 0.55901699437494742410;
  reldevx(x,y);

  x := qarg(a);
  y := 2.6318529747582863213;
  reldevx(x,y);

  x := qdot(a,b);
  y := -2.9850098891679916680;
  reldevx(x,y);

  qadd(a,b,c);     reldevq(c, raddab);
  qsub(a,b,c);     reldevq(c, rsubab);
  qmul(a,b,c);     reldevq(c, rmulab);
  qmul(b,a,c);     reldevq(c, rmulba);
  qdiv(a,b,c);     reldevq(c, rdivab);
  qconj(a,c);      reldevq(c, rconj);
  qneg(a,c);       reldevq(c, rneg);

  qarccos(a,c);    reldevq(c, racos);
  qarccosh(a,c);   reldevq(c, racosh);
  qarccot(a,c);    reldevq(c, racot);
  qarccoth(a,c);   reldevq(c, racoth);
  qarccsc(a,c);    reldevq(c, racsc);
  qarccsch(a,c);   reldevq(c, racsch);
  qarcsec(a,c);    reldevq(c, rasec);
  qarcsech(a,c);   reldevq(c, rasech);
  qarcsin(a,c);    reldevq(c, rasin);
  qarcsinh(a,c);   reldevq(c, rasinh);
  qarctan(a,c);    reldevq(c, ratan);
  qarctanh(a,c);   reldevq(c, ratanh);
  qcbrt(a,c);      reldevq(c, rcbrt);
  qcos(a,c);       reldevq(c, rcos);
  qcosh(a,c);      reldevq(c, rcosh);
  qcot(a,c);       reldevq(c, rcot);
  qcoth(a,c);      reldevq(c, rcoth);
  qcsc(a,c);       reldevq(c, rcsc);
  qcsch(a,c);      reldevq(c, rcsch);
  qexp(a,c);       reldevq(c, rexp);
  qinv(a,c);       reldevq(c, rinv);
  qln(a,c);        reldevq(c, rln);
  qpowx(a,5/2,c);  reldevq(c, rpowx);
  qsec(a,c);       reldevq(c, rsec);
  qsech(a,c);      reldevq(c, rsech);
  qsin(a,c);       reldevq(c, rsin);
  qsinh(a,c);      reldevq(c, rsinh);
  qsqrt(a,c);      reldevq(c, rsqrt);
  qtan(a,c);       reldevq(c, rtan);
  qtanh(a,c);      reldevq(c, rtanh);
  qunit(a,c);      reldevq(c, runit);
end;


{---------------------------------------------------------------------------}
function randx: double;
begin
  randx := 10.0*(random - 0.5);
end;


{---------------------------------------------------------------------------}
function reldev2q(const a,b: quaternion): double;
var
  t: quaternion;
  r: double;
begin
  inc(cnt);
  r := qabs(b);
  qsub(a,b,t);
  if r<0.5 then r := qabs(t)
  else r := qabs(t)/r;
  reldev2q := r;
  if r > eps_f then begin
    writeln(cnt:4, r);
    inc(failed);
  end;
  if r>maxerr_f then maxerr_f := r;
end;


{---------------------------------------------------------------------------}
procedure test_quat_func;
var
  i: integer;
begin
  for i:=1 to 100 do begin
    a.r := randx;
    a.x := randx;
    a.y := randx;
    a.z := randx;

    qln(a,b);
    qexp(b,c);
    reldev2q(c,a);

    qsqrt(a,b);
    qmul(b,b,c);
    reldev2q(c,a);

    qarccos(a,b);
    qcos(b,c);
    reldev2q(c,a);

    qarccosh(a,b);
    qcosh(b,c);
    reldev2q(c,a);

    qarccot(a,b);
    qcot(b,c);
    reldev2q(c,a);

    qarccoth(a,b);
    qcoth(b,c);
    reldev2q(c,a);

    qarccsc(a,b);
    qcsc(b,c);
    reldev2q(c,a);

    qarccsch(a,b);
    qcsch(b,c);
    reldev2q(c,a);

    qarcsec(a,b);
    qsec(b,c);
    reldev2q(c,a);

    qarcsech(a,b);
    qsech(b,c);
    reldev2q(c,a);

    qarcsin(a,b);
    qsin(b,c);
    reldev2q(c,a);

    qarcsinh(a,b);
    qsinh(b,c);
    reldev2q(c,a);

    qarctan(a,b);
    qtan(b,c);
    reldev2q(c,a);

    qarctanh(a,b);
    qtanh(b,c);
    reldev2q(c,a);
  end;
end;

{---------------------------------------------------------------------------}
procedure test_quat;
begin
  writeln;
  writeln('Test program for DAMQuat  (c) W.Ehrhardt 2017');
  eps := 12*eps_d;
  maxerr := -1;
  eps_f  := 32*eps_d;
  maxerr_f := -1;
  cnt := 0;
  failed := 0;
  test1;
  test2;
  test_quat_func;
  if failed=0 then writeln('All tests passed')
  else begin
    writeln(' *** number of failed tests: ', failed);
  end;
  writeln(' - max. error = ', maxerr/eps_d:1:2,',  max. error func(inv(q)) = ', maxerr_f/eps_d:1:2, ' [eps_d]');
end;

end.


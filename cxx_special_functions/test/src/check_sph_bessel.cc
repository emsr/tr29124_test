
//  sph_bessel

//  Compare against values generated by the GNU Scientific Library.
//  The GSL can be found on the web: http://www.gnu.org/software/gsl/
#include "verify.h"

// Test data for n=0.
// max(|f - f_GSL|): 3.3306690738754696e-16 at index 1
// max(|f - f_GSL| / |f_GSL|): 2.0843271082049370e-15
// mean(f - f_GSL): -3.0398963769498337e-17
// variance(f - f_GSL): 1.7755239868255660e-34
// stddev(f - f_GSL): 1.3324878936881814e-17
const testcase_sph_bessel<double>
data001[21] =
{
  { 1.0000000000000000, 0, 0.0000000000000000, 0.0 },
  { 0.98961583701809175, 0, 0.25000000000000000, 0.0 },
  { 0.95885107720840601, 0, 0.50000000000000000, 0.0 },
  { 0.90885168003111216, 0, 0.75000000000000000, 0.0 },
  { 0.84147098480789650, 0, 1.0000000000000000, 0.0 },
  { 0.75918769548446896, 0, 1.2500000000000000, 0.0 },
  { 0.66499665773603633, 0, 1.5000000000000000, 0.0 },
  { 0.56227768392796396, 0, 1.7500000000000000, 0.0 },
  { 0.45464871341284085, 0, 2.0000000000000000, 0.0 },
  { 0.34581030972796500, 0, 2.2500000000000000, 0.0 },
  { 0.23938885764158263, 0, 2.5000000000000000, 0.0 },
  { 0.13878581529175696, 0, 2.7500000000000000, 0.0 },
  { 0.047040002686622402, 0, 3.0000000000000000, 0.0 },
  { -0.033290810624648733, 0, 3.2500000000000000, 0.0 },
  { -0.10022377933989138, 0, 3.5000000000000000, 0.0 },
  { -0.15241635166462500, 0, 3.7500000000000000, 0.0 },
  { -0.18920062382698205, 0, 4.0000000000000000, 0.0 },
  { -0.21058573134790201, 0, 4.2500000000000000, 0.0 },
  { -0.21722891503668823, 0, 4.5000000000000000, 0.0 },
  { -0.21037742925797431, 0, 4.7500000000000000, 0.0 },
  { -0.19178485493262770, 0, 5.0000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for n=1.
// max(|f - f_GSL|): 3.1918911957973251e-16 at index 1
// max(|f - f_GSL| / |f_GSL|): 2.8516043985912409e-14
// mean(f - f_GSL): -1.9205867055457507e-17
// variance(f - f_GSL): 2.6407111150141003e-35
// stddev(f - f_GSL): 5.1387849877321193e-18
const testcase_sph_bessel<double>
data002[21] =
{
  { 0.0000000000000000, 1, 0.0000000000000000, 0.0 },
  { 0.082813661229788060, 1, 0.25000000000000000, 0.0 },
  { 0.16253703063606650, 1, 0.50000000000000000, 0.0 },
  { 0.23621708154305501, 1, 0.75000000000000000, 0.0 },
  { 0.30116867893975674, 1, 1.0000000000000000, 0.0 },
  { 0.35509226647136022, 1, 1.2500000000000000, 0.0 },
  { 0.39617297071222229, 1, 1.5000000000000000, 0.0 },
  { 0.42315642261568914, 1, 1.7500000000000000, 0.0 },
  { 0.43539777497999166, 1, 2.0000000000000000, 0.0 },
  { 0.43288174775586852, 1, 2.2500000000000000, 0.0 },
  { 0.41621298927540656, 1, 2.5000000000000000, 0.0 },
  { 0.38657752506335291, 1, 2.7500000000000000, 0.0 },
  { 0.34567749976235596, 1, 3.0000000000000000, 0.0 },
  { 0.29564272783258383, 1, 3.2500000000000000, 0.0 },
  { 0.23892368798597285, 1, 3.5000000000000000, 0.0 },
  { 0.17817146817998289, 1, 3.7500000000000000, 0.0 },
  { 0.11611074925915747, 1, 4.0000000000000000, 0.0 },
  { 0.055412178486091958, 1, 4.2500000000000000, 0.0 },
  { -0.0014295812457574522, 1, 4.5000000000000000, 0.0 },
  { -0.052206227820200179, 1, 4.7500000000000000, 0.0 },
  { -0.095089408079170795, 1, 5.0000000000000000, 0.0 },
};
const double toler002 = 2.5000000000000015e-12;

// Test data for n=2.
// max(|f - f_GSL|): 1.3877787807814457e-16 at index 9
// max(|f - f_GSL| / |f_GSL|): 8.2402809558653913e-16
// mean(f - f_GSL): -9.8300996972019074e-18
// variance(f - f_GSL): 5.3759274977605952e-34
// stddev(f - f_GSL): 2.3186046445568495e-17
const testcase_sph_bessel<double>
data003[21] =
{
  { 0.0000000000000000, 2, 0.0000000000000000, 0.0 },
  { 0.0041480977393611252, 2, 0.25000000000000000, 0.0 },
  { 0.016371106607993412, 2, 0.50000000000000000, 0.0 },
  { 0.036016646141108236, 2, 0.75000000000000000, 0.0 },
  { 0.062035052011373860, 2, 1.0000000000000000, 0.0 },
  { 0.093033744046795624, 2, 1.2500000000000000, 0.0 },
  { 0.12734928368840817, 2, 1.5000000000000000, 0.0 },
  { 0.16313332627036031, 2, 1.7500000000000000, 0.0 },
  { 0.19844794905714661, 2, 2.0000000000000000, 0.0 },
  { 0.23136535394652638, 2, 2.2500000000000000, 0.0 },
  { 0.26006672948890525, 2, 2.5000000000000000, 0.0 },
  { 0.28293512114099162, 2, 2.7500000000000000, 0.0 },
  { 0.29863749707573356, 2, 3.0000000000000000, 0.0 },
  { 0.30619179016241843, 2, 3.2500000000000000, 0.0 },
  { 0.30501551189929671, 2, 3.5000000000000000, 0.0 },
  { 0.29495352620861132, 2, 3.7500000000000000, 0.0 },
  { 0.27628368577135015, 2, 4.0000000000000000, 0.0 },
  { 0.24970021027926106, 2, 4.2500000000000000, 0.0 },
  { 0.21627586087284995, 2, 4.5000000000000000, 0.0 },
  { 0.17740507484521628, 2, 4.7500000000000000, 0.0 },
  { 0.13473121008512523, 2, 5.0000000000000000, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for n=5.
// max(|f - f_GSL|): 9.7144514654701197e-17 at index 18
// max(|f - f_GSL| / |f_GSL|): 2.7459190669103549e-15
// mean(f - f_GSL): 3.0672884393007950e-18
// variance(f - f_GSL): 4.9393356441808615e-37
// stddev(f - f_GSL): 7.0280407256794272e-19
const testcase_sph_bessel<double>
data004[21] =
{
  { 0.0000000000000000, 5, 0.0000000000000000, 0.0 },
  { 9.3719811237268220e-08, 5, 0.25000000000000000, 0.0 },
  { 2.9774668754574453e-06, 5, 0.50000000000000000, 0.0 },
  { 2.2339447678335762e-05, 5, 0.75000000000000000, 0.0 },
  { 9.2561158611258144e-05, 5, 1.0000000000000000, 0.0 },
  { 0.00027638888920123806, 5, 1.2500000000000000, 0.0 },
  { 0.00066962059628932456, 5, 1.5000000000000000, 0.0 },
  { 0.0014021729022572799, 5, 1.7500000000000000, 0.0 },
  { 0.0026351697702441169, 5, 2.0000000000000000, 0.0 },
  { 0.0045540034750567553, 5, 2.2500000000000000, 0.0 },
  { 0.0073576387377689359, 5, 2.5000000000000000, 0.0 },
  { 0.011244740276407145, 5, 2.7500000000000000, 0.0 },
  { 0.016397480955999105, 5, 3.0000000000000000, 0.0 },
  { 0.022964112474845508, 5, 3.2500000000000000, 0.0 },
  { 0.031041536537391189, 5, 3.5000000000000000, 0.0 },
  { 0.040659189440948935, 5, 3.7500000000000000, 0.0 },
  { 0.051765539757363456, 5, 4.0000000000000000, 0.0 },
  { 0.064218395773425613, 5, 4.2500000000000000, 0.0 },
  { 0.077780030832892866, 5, 4.5000000000000000, 0.0 },
  { 0.092117870593729223, 5, 4.7500000000000000, 0.0 },
  { 0.10681116145650453, 5, 5.0000000000000000, 0.0 },
};
const double toler004 = 2.5000000000000020e-13;

// Test data for n=10.
// max(|f - f_GSL|): 8.6736173798840355e-19 at index 18
// max(|f - f_GSL| / |f_GSL|): 6.7232224139500876e-15
// mean(f - f_GSL): 2.9233291056697913e-20
// variance(f - f_GSL): 4.2203838920690281e-38
// stddev(f - f_GSL): 2.0543572941601537e-19
const testcase_sph_bessel<double>
data005[21] =
{
  { 0.0000000000000000, 10, 0.0000000000000000, 0.0 },
  { 6.9267427453708468e-17, 10, 0.25000000000000000, 0.0 },
  { 7.0641239636618740e-14, 10, 0.50000000000000000, 0.0 },
  { 4.0459307474109287e-12, 10, 0.75000000000000000, 0.0 },
  { 7.1165526400473096e-11, 10, 1.0000000000000000, 0.0 },
  { 6.5470739530199939e-10, 10, 1.2500000000000000, 0.0 },
  { 3.9934406994836296e-09, 10, 1.5000000000000000, 0.0 },
  { 1.8327719460735247e-08, 10, 1.7500000000000000, 0.0 },
  { 6.8253008649747220e-08, 10, 2.0000000000000000, 0.0 },
  { 2.1653870546846626e-07, 10, 2.2500000000000000, 0.0 },
  { 6.0504362296385381e-07, 10, 2.5000000000000000, 0.0 },
  { 1.5246485352158441e-06, 10, 2.7500000000000000, 0.0 },
  { 3.5260038931752543e-06, 10, 3.0000000000000000, 0.0 },
  { 7.5839040020531456e-06, 10, 3.2500000000000000, 0.0 },
  { 1.5327786999397106e-05, 10, 3.5000000000000000, 0.0 },
  { 2.9348811002317661e-05, 10, 3.7500000000000000, 0.0 },
  { 5.3589865768632612e-05, 10, 4.0000000000000000, 0.0 },
  { 9.3818602410477989e-05, 10, 4.2500000000000000, 0.0 },
  { 0.00015817516371455801, 10, 4.5000000000000000, 0.0 },
  { 0.00025777607369970674, 10, 4.7500000000000000, 0.0 },
  { 0.00040734424424946052, 10, 5.0000000000000000, 0.0 },
};
const double toler005 = 5.0000000000000039e-13;

// Test data for n=20.
// max(|f - f_GSL|): 4.9275407583725281e-26 at index 20
// max(|f - f_GSL| / |f_GSL|): 2.4002866288153026e-14
// mean(f - f_GSL): 5.5808687000055343e-28
// variance(f - f_GSL): 1.2460231022000541e-52
// stddev(f - f_GSL): 1.1162540491304182e-26
const testcase_sph_bessel<double>
data006[21] =
{
  { 0.0000000000000000, 20, 0.0000000000000000, 0.0 },
  { 6.9307487073399339e-38, 20, 0.25000000000000000, 0.0 },
  { 7.2515880810153944e-32, 20, 0.50000000000000000, 0.0 },
  { 2.4025911398834722e-28, 20, 0.75000000000000000, 0.0 },
  { 7.5377957222368705e-26, 20, 1.0000000000000000, 0.0 },
  { 6.4953439243593413e-24, 20, 1.2500000000000000, 0.0 },
  { 2.4703120390884050e-22, 20, 1.5000000000000000, 0.0 },
  { 5.3404627138297197e-21, 20, 1.7500000000000000, 0.0 },
  { 7.6326411008876072e-20, 20, 2.0000000000000000, 0.0 },
  { 7.9496335952781075e-19, 20, 2.2500000000000000, 0.0 },
  { 6.4488532759578977e-18, 20, 2.5000000000000000, 0.0 },
  { 4.2725223040880135e-17, 20, 2.7500000000000000, 0.0 },
  { 2.3942249272752627e-16, 20, 3.0000000000000000, 0.0 },
  { 1.1654033741499860e-15, 20, 3.2500000000000000, 0.0 },
  { 5.0303402625237510e-15, 20, 3.5000000000000000, 0.0 },
  { 1.9572475798118559e-14, 20, 3.7500000000000000, 0.0 },
  { 6.9559880644906101e-14, 20, 4.0000000000000000, 0.0 },
  { 2.2825949745670935e-13, 20, 4.2500000000000000, 0.0 },
  { 6.9781823021792824e-13, 20, 4.5000000000000000, 0.0 },
  { 2.0024157388665026e-12, 20, 4.7500000000000000, 0.0 },
  { 5.4277267607932098e-12, 20, 5.0000000000000000, 0.0 },
};
const double toler006 = 2.5000000000000015e-12;
//  sph_bessel

// Test data for n=0.
// max(|f - f_GSL|): 1.0694570229397016e-15 at index 15
// max(|f - f_GSL| / |f_GSL|): 3.7496052611150890e-13
// mean(f - f_GSL): -2.6124109489412632e-17
// variance(f - f_GSL): 1.2558410443492295e-32
// stddev(f - f_GSL): 1.1206431387150995e-16
const testcase_sph_bessel<double>
data007[21] =
{
  { 1.0000000000000000, 0, 0.0000000000000000, 0.0 },
  { -0.19178485493262770, 0, 5.0000000000000000, 0.0 },
  { -0.054402111088936979, 0, 10.000000000000000, 0.0 },
  { 0.043352522677141118, 0, 15.000000000000000, 0.0 },
  { 0.045647262536381385, 0, 20.000000000000000, 0.0 },
  { -0.0052940700039109216, 0, 25.000000000000000, 0.0 },
  { -0.032934387469762058, 0, 30.000000000000000, 0.0 },
  { -0.012233790557032886, 0, 35.000000000000000, 0.0 },
  { 0.018627829011983722, 0, 40.000000000000000, 0.0 },
  { 0.018908967211869299, 0, 45.000000000000000, 0.0 },
  { -0.0052474970740785751, 0, 50.000000000000000, 0.0 },
  { -0.018177366788338544, 0, 55.000000000000000, 0.0 },
  { -0.0050801770183702783, 0, 60.000000000000000, 0.0 },
  { 0.012720441222924667, 0, 65.000000000000000, 0.0 },
  { 0.011055581165112701, 0, 70.000000000000000, 0.0 },
  { -0.0051704218054590724, 0, 75.000000000000000, 0.0 },
  { -0.012423608174042190, 0, 80.000000000000000, 0.0 },
  { -0.0020714778817480834, 0, 85.000000000000000, 0.0 },
  { 0.0099332962622284207, 0, 90.000000000000000, 0.0 },
  { 0.0071922285761696946, 0, 95.000000000000000, 0.0 },
  { -0.0050636564110975880, 0, 100.00000000000000, 0.0 },
};
const double toler007 = 2.5000000000000014e-11;

// Test data for n=1.
// max(|f - f_GSL|): 1.0044048925905713e-15 at index 13
// max(|f - f_GSL| / |f_GSL|): 6.5465850130521528e-13
// mean(f - f_GSL): 5.1561557602917773e-17
// variance(f - f_GSL): 2.9340234335545131e-33
// stddev(f - f_GSL): 5.4166626566129380e-17
const testcase_sph_bessel<double>
data008[21] =
{
  { 0.0000000000000000, 1, 0.0000000000000000, 0.0 },
  { -0.095089408079170795, 1, 5.0000000000000000, 0.0 },
  { 0.078466941798751549, 1, 10.000000000000000, 0.0 },
  { 0.053536029035730827, 1, 15.000000000000000, 0.0 },
  { -0.018121739963850528, 1, 20.000000000000000, 0.0 },
  { -0.039859875274695380, 1, 25.000000000000000, 0.0 },
  { -0.0062395279119115375, 1, 30.000000000000000, 0.0 },
  { 0.025470240415270681, 1, 35.000000000000000, 0.0 },
  { 0.017139147266606140, 1, 40.000000000000000, 0.0 },
  { -0.011253622702352454, 1, 45.000000000000000, 0.0 },
  { -0.019404270511323839, 1, 50.000000000000000, 0.0 },
  { -0.00073280223727807778, 1, 55.000000000000000, 0.0 },
  { 0.015788880056613101, 1, 60.000000000000000, 0.0 },
  { 0.0088488352686322581, 1, 65.000000000000000, 0.0 },
  { -0.0088894803131598157, 1, 70.000000000000000, 0.0 },
  { -0.012358955887069445, 1, 75.000000000000000, 0.0 },
  { 0.0012245454458125673, 1, 80.000000000000000, 0.0 },
  { 0.011556531358968161, 1, 85.000000000000000, 0.0 },
  { 0.0050889656932377614, 1, 90.000000000000000, 0.0 },
  { -0.0076103298149331573, 1, 95.000000000000000, 0.0 },
  { -0.0086738252869878150, 1, 100.00000000000000, 0.0 },
};
const double toler008 = 5.0000000000000028e-11;

// Test data for n=2.
// max(|f - f_GSL|): 1.0772632785815972e-15 at index 15
// max(|f - f_GSL| / |f_GSL|): 3.4761702917932150e-13
// mean(f - f_GSL): 2.4812741147453973e-17
// variance(f - f_GSL): 1.2315145072619553e-32
// stddev(f - f_GSL): 1.1097362331932554e-16
const testcase_sph_bessel<double>
data009[21] =
{
  { 0.0000000000000000, 2, 0.0000000000000000, 0.0 },
  { 0.13473121008512523, 2, 5.0000000000000000, 0.0 },
  { 0.077942193628562445, 2, 10.000000000000000, 0.0 },
  { -0.032645316869994959, 2, 15.000000000000000, 0.0 },
  { -0.048365523530958965, 2, 20.000000000000000, 0.0 },
  { 0.00051088497094747614, 2, 25.000000000000000, 0.0 },
  { 0.032310434678570907, 2, 30.000000000000000, 0.0 },
  { 0.014416954021198945, 2, 35.000000000000000, 0.0 },
  { -0.017342392966988262, 2, 40.000000000000000, 0.0 },
  { -0.019659208725359461, 2, 45.000000000000000, 0.0 },
  { 0.0040832408433991458, 2, 50.000000000000000, 0.0 },
  { 0.018137395757214285, 2, 55.000000000000000, 0.0 },
  { 0.0058696210212009327, 2, 60.000000000000000, 0.0 },
  { -0.012312033441295488, 2, 65.000000000000000, 0.0 },
  { -0.011436558892819550, 2, 70.000000000000000, 0.0 },
  { 0.0046760635699762939, 2, 75.000000000000000, 0.0 },
  { 0.012469528628260161, 2, 80.000000000000000, 0.0 },
  { 0.0024793554591234306, 2, 85.000000000000000, 0.0 },
  { -0.0097636640724538277, 2, 90.000000000000000, 0.0 },
  { -0.0074325547808517939, 2, 95.000000000000000, 0.0 },
  { 0.0048034416524879537, 2, 100.00000000000000, 0.0 },
};
const double toler009 = 2.5000000000000014e-11;

// Test data for n=5.
// max(|f - f_GSL|): 9.4455693266937146e-16 at index 14
// max(|f - f_GSL| / |f_GSL|): 8.4346477099300519e-13
// mean(f - f_GSL): 4.5743005943912237e-17
// variance(f - f_GSL): 1.4752048613244596e-33
// stddev(f - f_GSL): 3.8408395714016221e-17
const testcase_sph_bessel<double>
data010[21] =
{
  { 0.0000000000000000, 5, 0.0000000000000000, 0.0 },
  { 0.10681116145650453, 5, 5.0000000000000000, 0.0 },
  { -0.055534511621452155, 5, 10.000000000000000, 0.0 },
  { 0.065968007076521951, 5, 15.000000000000000, 0.0 },
  { 0.016683908063095682, 5, 20.000000000000000, 0.0 },
  { -0.036117795989722382, 5, 25.000000000000000, 0.0 },
  { -0.020504008736827489, 5, 30.000000000000000, 0.0 },
  { 0.018499481206814560, 5, 35.000000000000000, 0.0 },
  { 0.022448773791044995, 5, 40.000000000000000, 0.0 },
  { -0.0048552694845020138, 5, 45.000000000000000, 0.0 },
  { -0.020048300563664877, 5, 50.000000000000000, 0.0 },
  { -0.0052999924455565742, 5, 55.000000000000000, 0.0 },
  { 0.014151556281331407, 5, 60.000000000000000, 0.0 },
  { 0.011354588594416778, 5, 65.000000000000000, 0.0 },
  { -0.0064983781785323573, 5, 70.000000000000000, 0.0 },
  { -0.013089909320064257, 5, 75.000000000000000, 0.0 },
  { -0.00096200450071302446, 5, 80.000000000000000, 0.0 },
  { 0.011048668899130202, 5, 85.000000000000000, 0.0 },
  { 0.0065639581708136037, 5, 90.000000000000000, 0.0 },
  { -0.0064646119368202771, 5, 95.000000000000000, 0.0 },
  { -0.0092901489349075713, 5, 100.00000000000000, 0.0 },
};
const double toler010 = 5.0000000000000028e-11;

// Test data for n=10.
// max(|f - f_GSL|): 1.1999949645069563e-15 at index 13
// max(|f - f_GSL| / |f_GSL|): 2.9533832871668437e-12
// mean(f - f_GSL): 4.1910867550924778e-17
// variance(f - f_GSL): 1.5079763856266103e-32
// stddev(f - f_GSL): 1.2279968996811883e-16
const testcase_sph_bessel<double>
data011[21] =
{
  { 0.0000000000000000, 10, 0.0000000000000000, 0.0 },
  { 0.00040734424424946052, 10, 5.0000000000000000, 0.0 },
  { 0.064605154492564265, 10, 10.000000000000000, 0.0 },
  { 0.0018969790010883577, 10, 15.000000000000000, 0.0 },
  { 0.039686698644626366, 10, 20.000000000000000, 0.0 },
  { -0.036253285601128581, 10, 25.000000000000000, 0.0 },
  { -0.014529646403897799, 10, 30.000000000000000, 0.0 },
  { 0.026281264603993857, 10, 35.000000000000000, 0.0 },
  { 0.013124803182748323, 10, 40.000000000000000, 0.0 },
  { -0.017600831383728983, 10, 45.000000000000000, 0.0 },
  { -0.015039221463465955, 10, 50.000000000000000, 0.0 },
  { 0.0095256289349167390, 10, 55.000000000000000, 0.0 },
  { 0.015822719394008339, 10, 60.000000000000000, 0.0 },
  { -0.0019391391708249754, 10, 65.000000000000000, 0.0 },
  { -0.014293389028395012, 10, 70.000000000000000, 0.0 },
  { -0.0044210285031696227, 10, 75.000000000000000, 0.0 },
  { 0.010516146958338813, 10, 80.000000000000000, 0.0 },
  { 0.0086736275131325726, 10, 85.000000000000000, 0.0 },
  { -0.0052905066357239322, 10, 90.000000000000000, 0.0 },
  { -0.010258326955210768, 10, 95.000000000000000, 0.0 },
  { -0.00019565785971342414, 10, 100.00000000000000, 0.0 },
};
const double toler011 = 2.5000000000000017e-10;

// Test data for n=20.
// max(|f - f_GSL|): 8.5521867365656590e-16 at index 17
// max(|f - f_GSL| / |f_GSL|): 2.3231623379380350e-13
// mean(f - f_GSL): 7.4661920833706034e-18
// variance(f - f_GSL): 1.4597276928432959e-38
// stddev(f - f_GSL): 1.2081919106016625e-19
const testcase_sph_bessel<double>
data012[21] =
{
  { 0.0000000000000000, 20, 0.0000000000000000, 0.0 },
  { 5.4277267607932098e-12, 20, 5.0000000000000000, 0.0 },
  { 2.3083719613194670e-06, 20, 10.000000000000000, 0.0 },
  { 0.0015467058510412498, 20, 15.000000000000000, 0.0 },
  { 0.038324851639805160, 20, 20.000000000000000, 0.0 },
  { 0.028500071484154645, 20, 25.000000000000000, 0.0 },
  { -0.014711593353429081, 20, 30.000000000000000, 0.0 },
  { -0.010797653070264229, 20, 35.000000000000000, 0.0 },
  { 0.026535391837540293, 20, 40.000000000000000, 0.0 },
  { -0.011582959134716393, 20, 45.000000000000000, 0.0 },
  { -0.015785029898269291, 20, 50.000000000000000, 0.0 },
  { 0.013885519185862741, 20, 55.000000000000000, 0.0 },
  { 0.011112458964023273, 20, 60.000000000000000, 0.0 },
  { -0.011938384963927568, 20, 65.000000000000000, 0.0 },
  { -0.010117695207156904, 20, 70.000000000000000, 0.0 },
  { 0.0089871214102383232, 20, 75.000000000000000, 0.0 },
  { 0.010400578884991936, 20, 80.000000000000000, 0.0 },
  { -0.0055359020656326700, 20, 85.000000000000000, 0.0 },
  { -0.010639343320787521, 20, 90.000000000000000, 0.0 },
  { 0.0018051661455979529, 20, 95.000000000000000, 0.0 },
  { 0.010107671283873054, 20, 100.00000000000000, 0.0 },
};
const double toler012 = 2.5000000000000014e-11;

// Test data for n=50.
// max(|f - f_GSL|): 9.7377618121785581e-16 at index 17
// max(|f - f_GSL| / |f_GSL|): 2.0735742618499052e-12
// mean(f - f_GSL): 6.2784361339472195e-18
// variance(f - f_GSL): 5.0165160423422681e-33
// stddev(f - f_GSL): 7.0827367890825003e-17
const testcase_sph_bessel<double>
data013[21] =
{
  { 0.0000000000000000, 50, 0.0000000000000000, 0.0 },
  { 2.8574793504401511e-46, 50, 5.0000000000000000, 0.0 },
  { 2.2306960232186471e-31, 50, 10.000000000000000, 0.0 },
  { 7.6804716640080804e-23, 50, 15.000000000000000, 0.0 },
  { 5.6500807918725220e-17, 50, 20.000000000000000, 0.0 },
  { 1.2540416973758975e-12, 50, 25.000000000000000, 0.0 },
  { 2.6901637185735326e-09, 50, 30.000000000000000, 0.0 },
  { 1.0167148174422245e-06, 50, 35.000000000000000, 0.0 },
  { 9.3949174038179069e-05, 50, 40.000000000000000, 0.0 },
  { 0.0024888927213794561, 50, 45.000000000000000, 0.0 },
  { 0.018829107369282647, 50, 50.000000000000000, 0.0 },
  { 0.026373198438145489, 50, 55.000000000000000, 0.0 },
  { -0.021230978268739001, 50, 60.000000000000000, 0.0 },
  { 0.016539881802291313, 50, 65.000000000000000, 0.0 },
  { -0.015985416061436664, 50, 70.000000000000000, 0.0 },
  { 0.015462548984405590, 50, 75.000000000000000, 0.0 },
  { -0.010638570118081819, 50, 80.000000000000000, 0.0 },
  { 0.00046961239784540793, 50, 85.000000000000000, 0.0 },
  { 0.0096065882189920251, 50, 90.000000000000000, 0.0 },
  { -0.010613873910261154, 50, 95.000000000000000, 0.0 },
  { 0.00057971408822774927, 50, 100.00000000000000, 0.0 },
};
const double toler013 = 2.5000000000000017e-10;

// Test data for n=100.
// max(|f - f_GSL|): 3.1225022567582528e-17 at index 20
// max(|f - f_GSL| / |f_GSL|): 8.7701893132122237e-14
// mean(f - f_GSL): -1.2140001788067151e-18
// variance(f - f_GSL): 4.7284726903029159e-35
// stddev(f - f_GSL): 6.8763890889789794e-18
const testcase_sph_bessel<double>
data014[21] =
{
  { 0.0000000000000000, 100, 0.0000000000000000, 0.0 },
  { 5.5356503033889938e-120, 100, 5.0000000000000000, 0.0 },
  { 5.8320401820058771e-90, 100, 10.000000000000000, 0.0 },
  { 1.7406387750766626e-72, 100, 15.000000000000000, 0.0 },
  { 3.5152711125317012e-60, 100, 20.000000000000000, 0.0 },
  { 9.8455459353815965e-51, 100, 25.000000000000000, 0.0 },
  { 4.0888596744301583e-43, 100, 30.000000000000000, 0.0 },
  { 8.8975854911133939e-37, 100, 35.000000000000000, 0.0 },
  { 2.1513492547733828e-31, 100, 40.000000000000000, 0.0 },
  { 9.3673586994539108e-27, 100, 45.000000000000000, 0.0 },
  { 1.0190122629310471e-22, 100, 50.000000000000000, 0.0 },
  { 3.4887804977690388e-19, 100, 55.000000000000000, 0.0 },
  { 4.4442883425555593e-16, 100, 60.000000000000000, 0.0 },
  { 2.3832619568710723e-13, 100, 65.000000000000000, 0.0 },
  { 5.8948384175607987e-11, 100, 70.000000000000000, 0.0 },
  { 7.1884446357022277e-09, 100, 75.000000000000000, 0.0 },
  { 4.5247964400095002e-07, 100, 80.000000000000000, 0.0 },
  { 1.5096093228779032e-05, 100, 85.000000000000000, 0.0 },
  { 0.00026825172647807507, 100, 90.000000000000000, 0.0 },
  { 0.0024744308520581117, 100, 95.000000000000000, 0.0 },
  { 0.010880477011438350, 100, 100.00000000000000, 0.0 },
};
const double toler014 = 5.0000000000000029e-12;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_sph_bessel<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::sph_bessel(data[i].n, data[i].x);
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const Ret f0 = data[i].f0;
	    const Ret diff = f - f0;
	    if (std::abs(diff) > max_abs_diff)
	      max_abs_diff = std::abs(diff);
	    if (std::abs(f0) > Ret(10) * eps
	     && std::abs(f) > Ret(10) * eps)
	      {
		const Ret frac = diff / f0;
		if (std::abs(frac) > max_abs_frac)
		  max_abs_frac = std::abs(frac);
	      }
	  }
      }
    int num_errors = 0;
    VERIFY(max_abs_frac < toler);
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test(data001, toler001);
  num_errors += test(data002, toler002);
  num_errors += test(data003, toler003);
  num_errors += test(data004, toler004);
  num_errors += test(data005, toler005);
  num_errors += test(data006, toler006);
  num_errors += test(data007, toler007);
  num_errors += test(data008, toler008);
  num_errors += test(data009, toler009);
  num_errors += test(data010, toler010);
  num_errors += test(data011, toler011);
  num_errors += test(data012, toler012);
  num_errors += test(data013, toler013);
  num_errors += test(data014, toler014);
  return num_errors;
}
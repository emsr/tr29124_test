// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// { dg-skip-if "no extensions in strict dialects" { *-*-* } { "-std=c++*" } }

// Copyright (C) 2016-2018 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

//  fresnel_s

//  Compare against values generated by the GNU Scientific Library.
//  The GSL can be found on the web: http://www.gnu.org/software/gsl/
#include <limits>
#include <cmath>
#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << max_abs_frac \
	<< '\n'; \
    }
#else
#  include <testsuite_hooks.h>
#endif
#include <specfun_testcase.h>

// Test data.
// max(|f - f_GSL|): 1.7763568394002505e-15 at index 388
// max(|f - f_GSL| / |f_GSL|): 3.5144682605146659e-15
// mean(f - f_GSL): -5.4814551335435888e-18
// variance(f - f_GSL): 7.5492868464523123e-38
// stddev(f - f_GSL): 2.7475965581672125e-19
const testcase_fresnel_s<double>
data001[400] =
{
  { -0.48400633070510829, -19.899999999999999, 0.0 },
  { -0.48395467808739862, -19.800000000000001, 0.0 },
  { -0.48400149911370260, -19.699999999999999, 0.0 },
  { -0.48426660452027720, -19.600000000000001, 0.0 },
  { -0.48491377747941872, -19.500000000000000, 0.0 },
  { -0.48613909281472029, -19.399999999999999, 0.0 },
  { -0.48814634752103442, -19.300000000000001, 0.0 },
  { -0.49110465381924012, -19.199999999999999, 0.0 },
  { -0.49508580185270318, -19.100000000000001, 0.0 },
  { -0.49998522816706670, -19.000000000000000, 0.0 },
  { -0.50544113866605012, -18.899999999999999, 0.0 },
  { -0.51078069055407904, -18.800000000000001, 0.0 },
  { -0.51503609940649175, -18.699999999999999, 0.0 },
  { -0.51707863324524506, -18.600000000000001, 0.0 },
  { -0.51590229813703492, -18.500000000000000, 0.0 },
  { -0.51103958743392730, -18.399999999999999, 0.0 },
  { -0.50300681013324244, -18.300000000000001, 0.0 },
  { -0.49357730991366133, -18.199999999999999, 0.0 },
  { -0.48562173622205806, -18.100000000000001, 0.0 },
  { -0.48231616863711996, -18.000000000000000, 0.0 },
  { -0.48576890261541211, -17.899999999999999, 0.0 },
  { -0.49553539633619692, -17.800000000000001, 0.0 },
  { -0.50789526139161101, -17.699999999999999, 0.0 },
  { -0.51680884657947235, -17.600000000000001, 0.0 },
  { -0.51681175100068677, -17.500000000000000, 0.0 },
  { -0.50675220967147405, -17.399999999999999, 0.0 },
  { -0.49192297686929382, -17.300000000000001, 0.0 },
  { -0.48208003757999340, -17.199999999999999, 0.0 },
  { -0.48510206596727673, -17.100000000000001, 0.0 },
  { -0.49997937729695491, -17.000000000000000, 0.0 },
  { -0.51539764447709147, -16.899999999999999, 0.0 },
  { -0.51762428978876340, -16.800000000000001, 0.0 },
  { -0.50329846726239413, -16.699999999999999, 0.0 },
  { -0.48523935964401099, -16.600000000000001, 0.0 },
  { -0.48216841210237449, -16.500000000000000, 0.0 },
  { -0.49875837057440425, -16.399999999999999, 0.0 },
  { -0.51724734258551208, -16.300000000000001, 0.0 },
  { -0.51515475231295027, -16.199999999999999, 0.0 },
  { -0.49361889140909210, -16.100000000000001, 0.0 },
  { -0.48010572438090104, -16.000000000000000, 0.0 },
  { -0.49408939931225765, -15.899999999999999, 0.0 },
  { -0.51699614840547081, -15.800000000000001, 0.0 },
  { -0.51457780946000409, -15.699999999999999, 0.0 },
  { -0.48908932243831682, -15.600000000000000, 0.0 },
  { -0.48101678542100390, -15.500000000000000, 0.0 },
  { -0.50511339143415168, -15.399999999999999, 0.0 },
  { -0.52060088591333276, -15.300000000000001, 0.0 },
  { -0.49871387741832207, -15.199999999999999, 0.0 },
  { -0.47892213668560024, -15.100000000000000, 0.0 },
  { -0.49996997980970231, -15.000000000000000, 0.0 },
  { -0.52136079347254827, -14.899999999999999, 0.0 },
  { -0.49868073820357006, -14.800000000000001, 0.0 },
  { -0.47855793668137103, -14.699999999999999, 0.0 },
  { -0.50539037942758225, -14.600000000000000, 0.0 },
  { -0.52029395715602378, -14.500000000000000, 0.0 },
  { -0.48818436001775284, -14.399999999999999, 0.0 },
  { -0.48399092037687058, -14.300000000000001, 0.0 },
  { -0.51890751279314973, -14.199999999999999, 0.0 },
  { -0.50667250444726164, -14.100000000000000, 0.0 },
  { -0.47726375944182031, -14.000000000000000, 0.0 },
  { -0.50738195402708397, -13.899999999999999, 0.0 },
  { -0.51779703234472796, -13.800000000000001, 0.0 },
  { -0.47948494035222217, -13.699999999999999, 0.0 },
  { -0.49849019405657224, -13.600000000000000, 0.0 },
  { -0.52179926218777695, -13.500000000000000, 0.0 },
  { -0.48172388020281326, -13.399999999999999, 0.0 },
  { -0.49584282435016208, -13.300000000000001, 0.0 },
  { -0.52243698027430019, -13.199999999999999, 0.0 },
  { -0.48014633729836548, -13.100000000000000, 0.0 },
  { -0.49995388448191252, -13.000000000000000, 0.0 },
  { -0.51976048787627938, -12.899999999999999, 0.0 },
  { -0.47592559998979783, -12.800000000000001, 0.0 },
  { -0.51098198816812912, -12.699999999999999, 0.0 },
  { -0.50934679542387651, -12.600000000000000, 0.0 },
  { -0.47645404274109909, -12.500000000000000, 0.0 },
  { -0.52384763551119051, -12.399999999999999, 0.0 },
  { -0.48866392316693613, -12.300000000000001, 0.0 },
  { -0.49345748600389738, -12.199999999999999, 0.0 },
  { -0.52107102540974548, -12.100000000000000, 0.0 },
  { -0.47347456491993545, -12.000000000000000, 0.0 },
  { -0.52184956426284745, -11.900000000000000, 0.0 },
  { -0.49012717563678643, -11.799999999999999, 0.0 },
  { -0.49526026247089550, -11.699999999999999, 0.0 },
  { -0.51754095565565494, -11.600000000000000, 0.0 },
  { -0.47440277911222917, -11.500000000000000, 0.0 },
  { -0.52786202710467089, -11.400000000000000, 0.0 },
  { -0.47513851068918839, -11.299999999999999, 0.0 },
  { -0.51806001208333818, -11.199999999999999, 0.0 },
  { -0.49078143050743306, -11.100000000000000, 0.0 },
  { -0.49992388379375108, -11.000000000000000, 0.0 },
  { -0.50866137066057282, -10.900000000000000, 0.0 },
  { -0.48413995246468466, -10.799999999999999, 0.0 },
  { -0.52142029880840302, -10.699999999999999, 0.0 },
  { -0.47460051102087003, -10.600000000000000, 0.0 },
  { -0.52804040799812990, -10.500000000000000, 0.0 },
  { -0.47033321896953340, -10.399999999999999, 0.0 },
  { -0.53060780186343792, -10.299999999999999, 0.0 },
  { -0.46884960810453119, -10.199999999999999, 0.0 },
  { -0.53151256657876367, -10.100000000000000, 0.0 },
  { -0.46816997858488224, -10.000000000000000, 0.0 },
  { -0.53214917021179087, -9.8999999999999986, 0.0 },
  { -0.46757780186966436, -9.7999999999999989, 0.0 },
  { -0.53250259856426552, -9.6999999999999993, 0.0 },
  { -0.46785709076181747, -9.5999999999999996, 0.0 },
  { -0.53099984915139853, -9.5000000000000000, 0.0 },
  { -0.47134449224687508, -9.3999999999999986, 0.0 },
  { -0.52466585995281956, -9.2999999999999989, 0.0 },
  { -0.48135192933844551, -9.1999999999999993, 0.0 },
  { -0.51041329516968126, -9.0999999999999996, 0.0 },
  { -0.49986104562968459, -9.0000000000000000, 0.0 },
  { -0.48855154658234701, -8.8999999999999986, 0.0 },
  { -0.52294093275877707, -8.7999999999999989, 0.0 },
  { -0.46773905745666927, -8.6999999999999993, 0.0 },
  { -0.53692769034508359, -8.5999999999999996, 0.0 },
  { -0.46534124898107443, -8.5000000000000000, 0.0 },
  { -0.52428476977992999, -8.3999999999999986, 0.0 },
  { -0.49323233333085087, -8.2999999999999989, 0.0 },
  { -0.48588179985713581, -8.1999999999999993, 0.0 },
  { -0.53203939564156111, -8.0999999999999996, 0.0 },
  { -0.46021421439301446, -8.0000000000000000, 0.0 },
  { -0.53234203474777286, -7.8999999999999986, 0.0 },
  { -0.48964534042977748, -7.7999999999999989, 0.0 },
  { -0.48201416807598574, -7.6999999999999993, 0.0 },
  { -0.53885324338963214, -7.5999999999999996, 0.0 },
  { -0.46070123294683052, -7.5000000000000000, 0.0 },
  { -0.51606558037760442, -7.3999999999999986, 0.0 },
  { -0.51894732785814335, -7.2999999999999989, 0.0 },
  { -0.45725153058288015, -7.1999999999999993, 0.0 },
  { -0.53601735451077936, -7.0999999999999996, 0.0 },
  { -0.49970478945344637, -7.0000000000000000, 0.0 },
  { -0.46243950777859727, -6.8999999999999986, 0.0 },
  { -0.54363545682350578, -6.7999999999999989, 0.0 },
  { -0.49150144634508530, -6.6999999999999993, 0.0 },
  { -0.46306950114045575, -6.5999999999999996, 0.0 },
  { -0.54537645524323364, -6.5000000000000000, 0.0 },
  { -0.49649222154594042, -6.3999999999999986, 0.0 },
  { -0.45554543050439905, -6.2999999999999989, 0.0 },
  { -0.53982097881694135, -6.1999999999999993, 0.0 },
  { -0.51647708279510340, -6.0999999999999996, 0.0 },
  { -0.44696076123693029, -6.0000000000000000, 0.0 },
  { -0.51633069150415545, -5.8999999999999986, 0.0 },
  { -0.54604728378953227, -5.7999999999999989, 0.0 },
  { -0.45952838264767693, -5.6999999999999993, 0.0 },
  { -0.47003880651486152, -5.5999999999999996, 0.0 },
  { -0.55368406277902171, -5.5000000000000000, 0.0 },
  { -0.51403198870191358, -5.3999999999999986, 0.0 },
  { -0.44046778860409763, -5.2999999999999989, 0.0 },
  { -0.49687565586010518, -5.1999999999999993, 0.0 },
  { -0.56239007973300581, -5.0999999999999996, 0.0 },
  { -0.49919138191711693, -5.0000000000000000, 0.0 },
  { -0.43506736178749367, -4.8999999999999986, 0.0 },
  { -0.49675021895894833, -4.7999999999999989, 0.0 },
  { -0.56714546901226315, -4.6999999999999993, 0.0 },
  { -0.51619233694905442, -4.5999999999999996, 0.0 },
  { -0.43427297504870355, -4.5000000000000000, 0.0 },
  { -0.46226801641104603, -4.3999999999999986, 0.0 },
  { -0.55399588766657859, -4.2999999999999989, 0.0 },
  { -0.56319888839661081, -4.1999999999999993, 0.0 },
  { -0.47579825703282802, -4.0999999999999996, 0.0 },
  { -0.42051575424692844, -4.0000000000000000, 0.0 },
  { -0.47520240235068995, -3.8999999999999986, 0.0 },
  { -0.56561873979513244, -3.8000000000000007, 0.0 },
  { -0.57498034988747238, -3.6999999999999993, 0.0 },
  { -0.49230948911099942, -3.5999999999999979, 0.0 },
  { -0.41524801197243744, -3.5000000000000000, 0.0 },
  { -0.42964946444392793, -3.3999999999999986, 0.0 },
  { -0.51928608498206230, -3.3000000000000007, 0.0 },
  { -0.59334946461860383, -3.1999999999999993, 0.0 },
  { -0.58181586808587316, -3.0999999999999979, 0.0 },
  { -0.49631299896737496, -3.0000000000000000, 0.0 },
  { -0.41014058705671275, -2.8999999999999986, 0.0 },
  { -0.39152844354317162, -2.8000000000000007, 0.0 },
  { -0.45291748761671979, -2.6999999999999993, 0.0 },
  { -0.54998932315272153, -2.5999999999999979, 0.0 },
  { -0.61918175581959300, -2.5000000000000000, 0.0 },
  { -0.61968996494568307, -2.3999999999999986, 0.0 },
  { -0.55315164156070273, -2.3000000000000007, 0.0 },
  { -0.45570461212465635, -2.1999999999999993, 0.0 },
  { -0.37427335937810224, -2.0999999999999979, 0.0 },
  { -0.34341567836369830, -2.0000000000000000, 0.0 },
  { -0.37334731781698199, -1.8999999999999986, 0.0 },
  { -0.45093876926758242, -1.8000000000000007, 0.0 },
  { -0.54919594032156915, -1.6999999999999993, 0.0 },
  { -0.63888768350938285, -1.5999999999999979, 0.0 },
  { -0.69750496008209295, -1.5000000000000000, 0.0 },
  { -0.71352507736341197, -1.3999999999999986, 0.0 },
  { -0.68633328553465023, -1.3000000000000007, 0.0 },
  { -0.62340091854624902, -1.1999999999999993, 0.0 },
  { -0.53649791109681855, -1.0999999999999979, 0.0 },
  { -0.43825914739035476, -1.0000000000000000, 0.0 },
  { -0.33977634439313892, -0.89999999999999858, 0.0 },
  { -0.24934139305391545, -0.79999999999999716, 0.0 },
  { -0.17213645786347692, -0.69999999999999929, 0.0 },
  { -0.11054020735938577, -0.59999999999999787, 0.0 },
  { -0.064732432859999273, -0.50000000000000000, 0.0 },
  { -0.033359432660612816, -0.39999999999999858, 0.0 },
  { -0.014116998006576188, -0.29999999999999716, 0.0 },
  { -0.0041876091616567159, -0.19999999999999929, 0.0 },
  { -0.00052358954761217715, -0.099999999999997868, 0.0 },
  { 0.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.00052358954761223288, 0.10000000000000142, 0.0 },
  { 0.0041876091616569397, 0.20000000000000284, 0.0 },
  { 0.014116998006576686, 0.30000000000000071, 0.0 },
  { 0.033359432660613725, 0.40000000000000213, 0.0 },
  { 0.064732432859999273, 0.50000000000000000, 0.0 },
  { 0.11054020735938770, 0.60000000000000142, 0.0 },
  { 0.17213645786347945, 0.70000000000000284, 0.0 },
  { 0.24934139305391850, 0.80000000000000071, 0.0 },
  { 0.33977634439314230, 0.90000000000000213, 0.0 },
  { 0.43825914739035476, 1.0000000000000000, 0.0 },
  { 0.53649791109682166, 1.1000000000000014, 0.0 },
  { 0.62340091854625168, 1.2000000000000028, 0.0 },
  { 0.68633328553465023, 1.3000000000000007, 0.0 },
  { 0.71352507736341231, 1.4000000000000021, 0.0 },
  { 0.69750496008209295, 1.5000000000000000, 0.0 },
  { 0.63888768350937997, 1.6000000000000014, 0.0 },
  { 0.54919594032156549, 1.7000000000000028, 0.0 },
  { 0.45093876926758242, 1.8000000000000007, 0.0 },
  { 0.37334731781698011, 1.9000000000000021, 0.0 },
  { 0.34341567836369830, 2.0000000000000000, 0.0 },
  { 0.37427335937810458, 2.1000000000000014, 0.0 },
  { 0.45570461212465962, 2.2000000000000028, 0.0 },
  { 0.55315164156070273, 2.3000000000000007, 0.0 },
  { 0.61968996494568440, 2.4000000000000021, 0.0 },
  { 0.61918175581959300, 2.5000000000000000, 0.0 },
  { 0.54998932315271809, 2.6000000000000014, 0.0 },
  { 0.45291748761671680, 2.7000000000000028, 0.0 },
  { 0.39152844354317162, 2.8000000000000007, 0.0 },
  { 0.41014058705671497, 2.9000000000000021, 0.0 },
  { 0.49631299896737496, 3.0000000000000000, 0.0 },
  { 0.58181586808587515, 3.1000000000000014, 0.0 },
  { 0.59334946461860238, 3.2000000000000028, 0.0 },
  { 0.51928608498206230, 3.3000000000000007, 0.0 },
  { 0.42964946444392560, 3.4000000000000021, 0.0 },
  { 0.41524801197243744, 3.5000000000000000, 0.0 },
  { 0.49230948911100286, 3.6000000000000014, 0.0 },
  { 0.57498034988747404, 3.7000000000000028, 0.0 },
  { 0.56561873979513244, 3.8000000000000007, 0.0 },
  { 0.47520240235068667, 3.9000000000000021, 0.0 },
  { 0.42051575424692844, 4.0000000000000000, 0.0 },
  { 0.47579825703282963, 4.1000000000000014, 0.0 },
  { 0.56319888839661281, 4.2000000000000028, 0.0 },
  { 0.55399588766657748, 4.3000000000000007, 0.0 },
  { 0.46226801641104320, 4.4000000000000021, 0.0 },
  { 0.43427297504870355, 4.5000000000000000, 0.0 },
  { 0.51619233694905631, 4.6000000000000014, 0.0 },
  { 0.56714546901226259, 4.7000000000000028, 0.0 },
  { 0.49675021895894644, 4.8000000000000007, 0.0 },
  { 0.43506736178749372, 4.9000000000000021, 0.0 },
  { 0.49919138191711693, 5.0000000000000000, 0.0 },
  { 0.56239007973300581, 5.1000000000000014, 0.0 },
  { 0.49687565586010130, 5.2000000000000028, 0.0 },
  { 0.44046778860409791, 5.3000000000000007, 0.0 },
  { 0.51403198870191680, 5.4000000000000021, 0.0 },
  { 0.55368406277902171, 5.5000000000000000, 0.0 },
  { 0.47003880651485980, 5.6000000000000014, 0.0 },
  { 0.45952838264767909, 5.7000000000000028, 0.0 },
  { 0.54604728378953316, 5.8000000000000007, 0.0 },
  { 0.51633069150415178, 5.9000000000000021, 0.0 },
  { 0.44696076123693029, 6.0000000000000000, 0.0 },
  { 0.51647708279510518, 6.1000000000000014, 0.0 },
  { 0.53982097881693902, 6.2000000000000028, 0.0 },
  { 0.45554543050439805, 6.3000000000000007, 0.0 },
  { 0.49649222154594469, 6.4000000000000021, 0.0 },
  { 0.54537645524323364, 6.5000000000000000, 0.0 },
  { 0.46306950114045442, 6.6000000000000014, 0.0 },
  { 0.49150144634508858, 6.7000000000000028, 0.0 },
  { 0.54363545682350534, 6.8000000000000007, 0.0 },
  { 0.46243950777859505, 6.9000000000000021, 0.0 },
  { 0.49970478945344637, 7.0000000000000000, 0.0 },
  { 0.53601735451077859, 7.1000000000000014, 0.0 },
  { 0.45725153058287926, 7.2000000000000028, 0.0 },
  { 0.51894732785814501, 7.3000000000000007, 0.0 },
  { 0.51606558037760042, 7.4000000000000021, 0.0 },
  { 0.46070123294683052, 7.5000000000000000, 0.0 },
  { 0.53885324338963259, 7.6000000000000014, 0.0 },
  { 0.48201416807598257, 7.7000000000000028, 0.0 },
  { 0.48964534042977914, 7.8000000000000007, 0.0 },
  { 0.53234203474777086, 7.9000000000000021, 0.0 },
  { 0.46021421439301446, 8.0000000000000000, 0.0 },
  { 0.53203939564156244, 8.1000000000000014, 0.0 },
  { 0.48588179985713270, 8.2000000000000028, 0.0 },
  { 0.49323233333085303, 8.3000000000000007, 0.0 },
  { 0.52428476977992755, 8.4000000000000021, 0.0 },
  { 0.46534124898107443, 8.5000000000000000, 0.0 },
  { 0.53692769034508370, 8.6000000000000014, 0.0 },
  { 0.46773905745666733, 8.7000000000000028, 0.0 },
  { 0.52294093275877818, 8.8000000000000007, 0.0 },
  { 0.48855154658234412, 8.9000000000000021, 0.0 },
  { 0.49986104562968459, 9.0000000000000000, 0.0 },
  { 0.51041329516967926, 9.1000000000000014, 0.0 },
  { 0.48135192933844884, 9.2000000000000028, 0.0 },
  { 0.52466585995281889, 9.3000000000000007, 0.0 },
  { 0.47134449224687713, 9.4000000000000021, 0.0 },
  { 0.53099984915139853, 9.5000000000000000, 0.0 },
  { 0.46785709076181770, 9.6000000000000014, 0.0 },
  { 0.53250259856426496, 9.7000000000000028, 0.0 },
  { 0.46757780186966447, 9.8000000000000007, 0.0 },
  { 0.53214917021179087, 9.9000000000000021, 0.0 },
  { 0.46816997858488224, 10.000000000000000, 0.0 },
  { 0.53151256657876367, 10.100000000000001, 0.0 },
  { 0.46884960810453141, 10.200000000000003, 0.0 },
  { 0.53060780186343759, 10.300000000000001, 0.0 },
  { 0.47033321896953428, 10.400000000000002, 0.0 },
  { 0.52804040799812990, 10.500000000000000, 0.0 },
  { 0.47460051102087092, 10.600000000000001, 0.0 },
  { 0.52142029880840013, 10.700000000000003, 0.0 },
  { 0.48413995246468683, 10.800000000000001, 0.0 },
  { 0.50866137066057038, 10.900000000000002, 0.0 },
  { 0.49992388379375108, 11.000000000000000, 0.0 },
  { 0.49078143050743073, 11.100000000000001, 0.0 },
  { 0.51806001208334129, 11.200000000000003, 0.0 },
  { 0.47513851068918728, 11.300000000000001, 0.0 },
  { 0.52786202710467089, 11.400000000000002, 0.0 },
  { 0.47440277911222917, 11.500000000000000, 0.0 },
  { 0.51754095565565439, 11.600000000000001, 0.0 },
  { 0.49526026247089927, 11.700000000000003, 0.0 },
  { 0.49012717563678504, 11.800000000000001, 0.0 },
  { 0.52184956426284790, 11.900000000000002, 0.0 },
  { 0.47347456491993545, 12.000000000000000, 0.0 },
  { 0.52107102540974448, 12.100000000000001, 0.0 },
  { 0.49345748600390099, 12.200000000000003, 0.0 },
  { 0.48866392316693280, 12.300000000000004, 0.0 },
  { 0.52384763551119051, 12.399999999999999, 0.0 },
  { 0.47645404274109909, 12.500000000000000, 0.0 },
  { 0.50934679542387451, 12.600000000000001, 0.0 },
  { 0.51098198816813234, 12.700000000000003, 0.0 },
  { 0.47592559998979678, 12.800000000000004, 0.0 },
  { 0.51976048787627938, 12.899999999999999, 0.0 },
  { 0.49995388448191252, 13.000000000000000, 0.0 },
  { 0.48014633729836470, 13.100000000000001, 0.0 },
  { 0.52243698027429863, 13.200000000000003, 0.0 },
  { 0.49584282435016613, 13.300000000000004, 0.0 },
  { 0.48172388020281326, 13.399999999999999, 0.0 },
  { 0.52179926218777695, 13.500000000000000, 0.0 },
  { 0.49849019405657491, 13.600000000000001, 0.0 },
  { 0.47948494035222033, 13.700000000000003, 0.0 },
  { 0.51779703234472541, 13.800000000000004, 0.0 },
  { 0.50738195402708397, 13.899999999999999, 0.0 },
  { 0.47726375944182031, 14.000000000000000, 0.0 },
  { 0.50667250444726042, 14.100000000000001, 0.0 },
  { 0.51890751279315117, 14.200000000000003, 0.0 },
  { 0.48399092037687325, 14.300000000000004, 0.0 },
  { 0.48818436001775284, 14.399999999999999, 0.0 },
  { 0.52029395715602378, 14.500000000000000, 0.0 },
  { 0.50539037942758347, 14.600000000000001, 0.0 },
  { 0.47855793668137159, 14.700000000000003, 0.0 },
  { 0.49868073820356640, 14.800000000000004, 0.0 },
  { 0.52136079347254827, 14.899999999999999, 0.0 },
  { 0.49996997980970231, 15.000000000000000, 0.0 },
  { 0.47892213668560024, 15.100000000000001, 0.0 },
  { 0.49871387741831974, 15.200000000000003, 0.0 },
  { 0.52060088591333231, 15.300000000000004, 0.0 },
  { 0.50511339143415168, 15.399999999999999, 0.0 },
  { 0.48101678542100390, 15.500000000000000, 0.0 },
  { 0.48908932243831482, 15.600000000000001, 0.0 },
  { 0.51457780946000087, 15.700000000000003, 0.0 },
  { 0.51699614840547270, 15.800000000000004, 0.0 },
  { 0.49408939931225765, 15.899999999999999, 0.0 },
  { 0.48010572438090104, 16.000000000000000, 0.0 },
  { 0.49361889140909210, 16.100000000000001, 0.0 },
  { 0.51515475231294738, 16.200000000000003, 0.0 },
  { 0.51724734258551364, 16.300000000000004, 0.0 },
  { 0.49875837057440425, 16.399999999999999, 0.0 },
  { 0.48216841210237449, 16.500000000000000, 0.0 },
  { 0.48523935964401099, 16.600000000000001, 0.0 },
  { 0.50329846726239091, 16.700000000000003, 0.0 },
  { 0.51762428978876229, 16.800000000000004, 0.0 },
  { 0.51539764447709147, 16.899999999999999, 0.0 },
  { 0.49997937729695491, 17.000000000000000, 0.0 },
  { 0.48510206596727673, 17.100000000000001, 0.0 },
  { 0.48208003757999263, 17.200000000000003, 0.0 },
  { 0.49192297686929098, 17.300000000000004, 0.0 },
  { 0.50675220967147405, 17.399999999999999, 0.0 },
  { 0.51681175100068677, 17.500000000000000, 0.0 },
  { 0.51680884657947235, 17.600000000000001, 0.0 },
  { 0.50789526139161367, 17.700000000000003, 0.0 },
  { 0.49553539633619986, 17.800000000000004, 0.0 },
  { 0.48576890261541211, 17.899999999999999, 0.0 },
  { 0.48231616863711996, 18.000000000000000, 0.0 },
  { 0.48562173622205806, 18.100000000000001, 0.0 },
  { 0.49357730991365767, 18.200000000000003, 0.0 },
  { 0.50300681013323856, 18.300000000000004, 0.0 },
  { 0.51103958743392131, 18.400000000000006, 0.0 },
  { 0.51590229813703492, 18.500000000000000, 0.0 },
  { 0.51707863324524506, 18.600000000000001, 0.0 },
  { 0.51503609940649264, 18.700000000000003, 0.0 },
  { 0.51078069055408204, 18.800000000000004, 0.0 },
  { 0.50544113866605733, 18.900000000000006, 0.0 },
  { 0.49998522816706670, 19.000000000000000, 0.0 },
  { 0.49508580185270318, 19.100000000000001, 0.0 },
  { 0.49110465381924173, 19.200000000000003, 0.0 },
  { 0.48814634752103703, 19.300000000000004, 0.0 },
  { 0.48613909281472428, 19.400000000000006, 0.0 },
  { 0.48491377747941872, 19.500000000000000, 0.0 },
  { 0.48426660452027720, 19.600000000000001, 0.0 },
  { 0.48400149911370283, 19.700000000000003, 0.0 },
  { 0.48395467808739889, 19.800000000000004, 0.0 },
  { 0.48400633070510840, 19.900000000000006, 0.0 },
  { 0.48408453592595391, 20.000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_fresnel_s<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::fresnel_s(data[i].x);
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
    VERIFY(max_abs_frac < toler);
  }

int
main()
{
  test(data001, toler001);
  return 0;
}
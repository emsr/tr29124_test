
//  ldouble_factorial

#include "verify.h"

// Test data.
// max(|f - f_GSL|): 1.9850858734571375e-09 at index 368
// max(|f - f_GSL| / |f_GSL|): 2.3061783919300212e-12
// mean(f - f_GSL): 1.9880565944871425e-10
// variance(f - f_GSL): 3.3134948081371490e-22
// stddev(f - f_GSL): 1.8203007466177528e-11
const testcase_ldouble_factorial<double>
data001[501] =
{
  { 0.0000000000000000, 0, 0.0 },
  { 0.0000000000000000, 1, 0.0 },
  { 0.69314718055994529, 2, 0.0 },
  { 1.0986122886681098, 3, 0.0 },
  { 2.0794415416798357, 4, 0.0 },
  { 2.7080502011022101, 5, 0.0 },
  { 3.8712010109078911, 6, 0.0 },
  { 4.6539603501575231, 7, 0.0 },
  { 5.9506425525877269, 8, 0.0 },
  { 6.8511849274937431, 9, 0.0 },
  { 8.2532276455817719, 10, 0.0 },
  { 9.2490802002921129, 11, 0.0 },
  { 10.738134295369774, 12, 0.0 },
  { 11.814029557753649, 13, 0.0 },
  { 13.377191624985031, 14, 0.0 },
  { 14.522079758855860, 15, 0.0 },
  { 16.149780347224812, 16, 0.0 },
  { 17.355293102912075, 17, 0.0 },
  { 19.040152105120978, 18, 0.0 },
  { 20.299732082078517, 19, 0.0 },
  { 22.035884378674968, 20, 0.0 },
  { 23.344254519801940, 21, 0.0 },
  { 25.126926832033284, 22, 0.0 },
  { 26.479748735731089, 23, 0.0 },
  { 28.304980662381229, 24, 0.0 },
  { 29.698624560599288, 25, 0.0 },
  { 31.563077200402713, 26, 0.0 },
  { 32.994461426603621, 27, 0.0 },
  { 34.895281710577919, 28, 0.0 },
  { 36.361757256590096, 29, 0.0 },
  { 38.296479092240070, 30, 0.0 },
  { 39.795744461075238, 31, 0.0 },
  { 41.762214995039798, 32, 0.0 },
  { 43.292252022541717, 33, 0.0 },
  { 45.288575519655957, 34, 0.0 },
  { 46.847600084031136, 35, 0.0 },
  { 48.872094458112066, 36, 0.0 },
  { 50.458517996675354, 37, 0.0 },
  { 52.509680617838455, 38, 0.0 },
  { 54.122079642805005, 39, 0.0 },
  { 56.198560071952393, 40, 0.0 },
  { 57.835651709509314, 41, 0.0 },
  { 59.936229690235763, 42, 0.0 },
  { 61.596851825202876, 43, 0.0 },
  { 63.720419324154022, 44, 0.0 },
  { 65.403514314973194, 45, 0.0 },
  { 67.549060720643112, 46, 0.0 },
  { 69.253661916683257, 47, 0.0 },
  { 71.420261731551008, 48, 0.0 },
  { 73.145482214793873, 49, 0.0 },
  { 75.332284736979148, 50, 0.0 },
  { 77.077307847518199, 51, 0.0 },
  { 79.283528455560585, 52, 0.0 },
  { 81.047599761070330, 53, 0.0 },
  { 83.272512502124854, 54, 0.0 },
  { 85.054932946302799, 55, 0.0 },
  { 87.297864192860004, 56, 0.0 },
  { 89.097984214137341, 57, 0.0 },
  { 91.358307203406426, 58, 0.0 },
  { 93.175521658043067, 59, 0.0 },
  { 95.452651765628531, 60, 0.0 },
  { 97.286395522216381, 61, 0.0 },
  { 99.579786150673613, 62, 0.0 },
  { 101.42953024860792, 63, 0.0 },
  { 103.73866923403328, 64, 0.0 },
  { 105.60391751850355, 65, 0.0 },
  { 107.92832397605972, 66, 0.0 },
  { 109.80861013789452, 67, 0.0 },
  { 112.14783168123581, 68, 0.0 },
  { 114.04271664249177, 69, 0.0 },
  { 116.39632692328517, 70, 0.0 },
  { 118.30539651953309, 71, 0.0 },
  { 120.67299304230123, 72, 0.0 },
  { 122.59585596068148, 73, 0.0 },
  { 124.97705813550540, 74, 0.0 },
  { 126.91334407421779, 75, 0.0 },
  { 129.30779147579173, 76, 0.0 },
  { 131.25714949607146, 77, 0.0 },
  { 133.66450030248131, 78, 0.0 },
  { 135.62659734853850, 79, 0.0 },
  { 138.04652693715522, 80, 0.0 },
  { 140.02104650321093, 81, 0.0 },
  { 142.45324618441947, 82, 0.0 },
  { 144.43988711100752, 83, 0.0 },
  { 146.88406298326277, 84, 0.0 },
  { 148.88253836749786, 85, 0.0 },
  { 151.33841027951630, 86, 0.0 },
  { 153.34844648615243, 87, 0.0 },
  { 155.81574709399450, 88, 0.0 },
  { 157.83708285588457, 89, 0.0 },
  { 160.31555676432475, 90, 0.0 },
  { 162.34794236240143, 91, 0.0 },
  { 164.83734534137380, 92, 0.0 },
  { 166.88054185555467, 93, 0.0 },
  { 169.38064012364379, 94, 0.0 },
  { 171.43441874715521, 95, 0.0 },
  { 173.94498831511163, 96, 0.0 },
  { 176.00912972565860, 97, 0.0 },
  { 178.52995579378219, 98, 0.0 },
  { 180.60424957579320, 99, 0.0 },
  { 183.13512597977029, 100, 0.0 },
  { 185.21937009263445, 101, 0.0 },
  { 187.76009879305457, 102, 0.0 },
  { 189.85409908086407, 103, 0.0 },
  { 192.40448969219594, 104, 0.0 },
  { 194.50805943102162, 105, 0.0 },
  { 197.06792878630802, 106, 0.0 },
  { 199.18088826548350, 107, 0.0 },
  { 201.75006001343223, 108, 0.0 },
  { 203.87223614771267, 109, 0.0 },
  { 206.45054037922463, 110, 0.0 },
  { 208.58176634902500, 111, 0.0 },
  { 211.16903925051975, 112, 0.0 },
  { 213.30915416773735, 113, 0.0 },
  { 215.90523769891422, 114, 0.0 },
  { 218.05408629610059, 115, 0.0 },
  { 220.65882789002060, 116, 0.0 },
  { 222.81626023089834, 117, 0.0 },
  { 225.42951251448628, 118, 0.0 },
  { 227.59538372400988, 119, 0.0 },
  { 230.21700425726831, 120, 0.0 },
  { 232.39117426960661, 121, 0.0 },
  { 235.02102530200156, 122, 0.0 },
  { 237.20335862497902, 123, 0.0 },
  { 239.84130686760659, 124, 0.0 },
  { 242.03167236228134, 125, 0.0 },
  { 244.67758877455807, 126, 0.0 },
  { 246.87585944873993, 127, 0.0 },
  { 249.52961903847770, 128, 0.0 },
  { 251.73567185310159, 129, 0.0 },
  { 254.39715348893327, 130, 0.0 },
  { 256.61086917630274, 131, 0.0 },
  { 259.27995541151967, 132, 0.0 },
  { 261.50121830452451, 133, 0.0 },
  { 264.17779521147054, 134, 0.0 },
  { 266.40649308296292, 135, 0.0 },
  { 269.09045009720660, 136, 0.0 },
  { 271.32647400879108, 137, 0.0 },
  { 274.01770378236381, 138, 0.0 },
  { 276.26094794192176, 139, 0.0 },
  { 278.95934620497314, 140, 0.0 },
  { 281.20970783229990, 141, 0.0 },
  { 283.91517326257440, 142, 0.0 },
  { 286.17255246255985, 143, 0.0 },
  { 288.88498656215040, 144, 0.0 },
  { 291.14928620498040, 145, 0.0 },
  { 293.86859318385871, 146, 0.0 },
  { 296.13971879175915, 147, 0.0 },
  { 298.86580545762286, 148, 0.0 },
  { 301.14366509770457, 149, 0.0 },
  { 303.87644075171909, 150, 0.0 },
  { 306.16094493451953, 151, 0.0 },
  { 308.90032127256535, 152, 0.0 },
  { 311.19138285591197, 153, 0.0 },
  { 313.93727387497898, 154, 0.0 },
  { 316.23480797283122, 155, 0.0 },
  { 318.98712988222854, 156, 0.0 },
  { 321.29105377817950, 157, 0.0 },
  { 324.04972491525552, 158, 0.0 },
  { 326.35995798039971, 159, 0.0 },
  { 329.12489873048935, 160, 0.0 },
  { 331.44136234538422, 161, 0.0 },
  { 334.21249506572173, 162, 0.0 },
  { 336.53511254619099, 163, 0.0 },
  { 339.31236149354589, 164, 0.0 },
  { 341.64105802009152, 165, 0.0 },
  { 344.42434928190244, 166, 0.0 },
  { 346.75905183250831, 167, 0.0 },
  { 349.54831326130574, 168, 0.0 },
  { 351.88895054743136, 169, 0.0 },
  { 354.68411169835599, 170, 0.0 },
  { 357.03061410393406, 171, 0.0 },
  { 359.83160617516944, 172, 0.0 },
  { 362.18390569843183, 173, 0.0 },
  { 364.99066147438396, 174, 0.0 },
  { 367.34869167235530, 175, 0.0 },
  { 370.16114546942208, 176, 0.0 },
  { 372.52484140492913, 177, 0.0 },
  { 375.34292901971418, 178, 0.0 },
  { 377.71222721076992, 179, 0.0 },
  { 380.53588587060443, 180, 0.0 },
  { 382.91072424203571, 181, 0.0 },
  { 385.73989255768117, 182, 0.0 },
  { 388.12021039487718, 183, 0.0 },
  { 390.95482831529017, 184, 0.0 },
  { 393.34056621995546, 185, 0.0 },
  { 396.18057498900339, 186, 0.0 },
  { 398.57167483681008, 187, 0.0 },
  { 401.41701695183332, 188, 0.0 },
  { 403.81342185186969, 189, 0.0 },
  { 406.66404102399383, 190, 0.0 },
  { 409.06569527991633, 191, 0.0 },
  { 411.92153639602162, 192, 0.0 },
  { 414.32838546882124, 193, 0.0 },
  { 417.18939455508496, 194, 0.0 },
  { 419.60138502738499, 195, 0.0 },
  { 422.46750921431544, 196, 0.0 },
  { 424.88458875612298, 197, 0.0 },
  { 427.75577624501000, 198, 0.0 },
  { 430.17789358084747, 199, 0.0 },
  { 433.05409361155802, 200, 0.0 },
  { 435.48119848890650, 201, 0.0 },
  { 438.36236130895924, 202, 0.0 },
  { 440.79440446794831, 203, 0.0 },
  { 443.68048130280346, 204, 0.0 },
  { 446.11741444708673, 205, 0.0 },
  { 449.00835747159300, 206, 0.0 },
  { 451.45013324035210, 207, 0.0 },
  { 454.34589555129435, 208, 0.0 },
  { 456.79246749231692, 209, 0.0 },
  { 459.69300308201178, 210, 0.0 },
  { 462.14432562579299, 211, 0.0 },
  { 465.04958935668384, 212, 0.0 },
  { 467.50561779150240, 213, 0.0 },
  { 470.41556537170567, 214, 0.0 },
  { 472.87625581963005, 215, 0.0 },
  { 475.79084377938983, 216, 0.0 },
  { 478.25615317317050, 217, 0.0 },
  { 481.17533884217892, 218, 0.0 },
  { 483.64522490298702, 219, 0.0 },
  { 486.56896638853129, 220, 0.0 },
  { 489.04338760450474, 221, 0.0 },
  { 491.97164377040355, 222, 0.0 },
  { 494.45055937596487, 223, 0.0 },
  { 497.38328982225863, 224, 0.0 },
  { 499.86665977816932, 225, 0.0 },
  { 502.80382482153090, 226, 0.0 },
  { 505.29160979565069, 227, 0.0 },
  { 508.23317045048532, 228, 0.0 },
  { 510.72533179920498, 229, 0.0 },
  { 513.67124975940851, 230, 0.0 },
  { 516.16774950972672, 231, 0.0 },
  { 519.11798713107487, 232, 0.0 },
  { 521.61878796329245, 233, 0.0 },
  { 524.57330824643259, 234, 0.0 },
  { 527.07837347743657, 235, 0.0 },
  { 530.03714005145810, 236, 0.0 },
  { 532.54643361857177, 237, 0.0 },
  { 535.50941072512967, 238, 0.0 },
  { 538.02289717050326, 239, 0.0 },
  { 540.99004964847165, 240, 0.0 },
  { 543.50769410399391, 241, 0.0 },
  { 546.47898737462833, 242, 0.0 },
  { 549.00075554733451, 243, 0.0 },
  { 551.97615559992153, 244, 0.0 },
  { 554.50201375787913, 245, 0.0 },
  { 557.48148713585385, 246, 0.0 },
  { 560.01140209450716, 247, 0.0 },
  { 562.99491588201886, 248, 0.0 },
  { 565.52885499097181, 249, 0.0 },
  { 568.51637679988107, 250, 0.0 },
  { 571.05430793010362, 251, 0.0 },
  { 574.04580588739248, 252, 0.0 },
  { 576.58769741883111, 253, 0.0 },
  { 579.58314015441101, 254, 0.0 },
  { 582.12896096398958, 255, 0.0 },
  { 585.12831759889059, 256, 0.0 },
  { 587.67803704888479, 257, 0.0 },
  { 590.68127718381220, 258, 0.0 },
  { 593.23486511058434, 259, 0.0 },
  { 596.24195881482774, 260, 0.0 },
  { 598.79938551790701, 261, 0.0 },
  { 601.81030331858881, 262, 0.0 },
  { 604.37153955008478, 263, 0.0 },
  { 607.38625242173521, 264, 0.0 },
  { 609.95126937607108, 265, 0.0 },
  { 612.96974873051693, 266, 0.0 },
  { 615.53851803447128, 267, 0.0 },
  { 618.56073571102775, 268, 0.0 },
  { 621.13322941407307, 269, 0.0 },
  { 624.15915767002616, 270, 0.0 },
  { 626.73534823495277, 271, 0.0 },
  { 629.76495973632211, 272, 0.0 },
  { 632.34482003013773, 273, 0.0 },
  { 635.37808784271022, 274, 0.0 },
  { 637.96159112780435, 275, 0.0 },
  { 640.99848870842732, 276, 0.0 },
  { 643.58560863399168, 277, 0.0 },
  { 646.62610982211800, 278, 0.0 },
  { 649.21682041581300, 279, 0.0 },
  { 652.26089942528722, 280, 0.0 },
  { 654.85517508514681, 281, 0.0 },
  { 657.90280649622537, 282, 0.0 },
  { 660.50062198278999, 283, 0.0 },
  { 663.55178073438651, 284, 0.0 },
  { 666.15311116305872, 285, 0.0 },
  { 669.20777254520635, 286, 0.0 },
  { 671.81259337881829, 287, 0.0 },
  { 674.87073302534236, 288, 0.0 },
  { 677.47902006693073, 289, 0.0 },
  { 680.54061394832286, 290, 0.0 },
  { 683.15234333410228, 291, 0.0 },
  { 686.21736775059117, 292, 0.0 },
  { 688.83251594311935, 293, 0.0 },
  { 691.90094751792981, 294, 0.0 },
  { 694.51949129945910, 295, 0.0 },
  { 697.59130697225385, 296, 0.0 },
  { 700.21322343826182, 297, 0.0 },
  { 703.28840045875927, 298, 0.0 },
  { 705.91366701165236, 299, 0.0 },
  { 708.99218293341528, 300, 0.0 },
  { 711.62077727640133, 301, 0.0 },
  { 714.70260995079025, 302, 0.0 },
  { 717.33451008191059, 303, 0.0 },
  { 720.41963765219657, 304, 0.0 },
  { 723.05482185851804, 305, 0.0 },
  { 726.14322275414884, 306, 0.0 },
  { 728.78166960610531, 307, 0.0 },
  { 731.87332253712248, 308, 0.0 },
  { 734.51501088300301, 309, 0.0 },
  { 737.60989483460150, 310, 0.0 },
  { 740.25480379518228, 311, 0.0 },
  { 743.35289802241107, 312, 0.0 },
  { 746.00100698572237, 313, 0.0 },
  { 749.10229100831930, 314, 0.0 },
  { 751.75357962454791, 315, 0.0 },
  { 754.85803322190623, 316, 0.0 },
  { 757.51248139842528, 317, 0.0 },
  { 760.62008460468644, 318, 0.0 },
  { 763.27767250121019, 319, 0.0 },
  { 766.38840560048016, 320, 0.0 },
  { 769.04911362434018, 321, 0.0 },
  { 772.16295714602461, 322, 0.0 },
  { 774.82676594756288, 323, 0.0 },
  { 777.94370066181682, 324, 0.0 },
  { 780.61059112989255, 325, 0.0 },
  { 783.73059804318370, 326, 0.0 },
  { 786.40055130078997, 327, 0.0 },
  { 789.52361165156788, 328, 0.0 },
  { 792.19660905155513, 329, 0.0 },
  { 795.32270430602841, 330, 0.0 },
  { 797.99872742693231, 331, 0.0 },
  { 801.12783927494479, 332, 0.0 },
  { 803.80686991691266, 333, 0.0 },
  { 806.93898026792147, 334, 0.0 },
  { 809.62100044873785, 335, 0.0 },
  { 812.75609142788471, 336, 0.0 },
  { 815.44108337909017, 337, 0.0 },
  { 818.57913732336772, 338, 0.0 },
  { 821.26708348647048, 339, 0.0 },
  { 824.40808294097792, 340, 0.0 },
  { 827.09896596375404, 341, 0.0 },
  { 830.24289367804067, 342, 0.0 },
  { 832.93669641092004, 343, 0.0 },
  { 836.08353533541401, 344, 0.0 },
  { 838.78024082795127, 345, 0.0 },
  { 841.92997411047179, 346, 0.0 },
  { 844.62956560789826, 347, 0.0 },
  { 847.78217659024608, 348, 0.0 },
  { 850.48463753010071, 349, 0.0 },
  { 853.64010974472967, 350, 0.0 },
  { 856.34542375356659, 351, 0.0 },
  { 859.50374092032780, 352, 0.0 },
  { 862.21189181049988, 353, 0.0 },
  { 865.37303783346147, 354, 0.0 },
  { 868.08400959997539, 355, 0.0 },
  { 871.24796856431340, 356, 0.0 },
  { 873.96174538175489, 357, 0.0 },
  { 877.12850155071419, 358, 0.0 },
  { 879.84506777024319, 359, 0.0 },
  { 883.01460558216434, 360, 0.0 },
  { 885.73394572857603, 361, 0.0 },
  { 888.90624979399013, 362, 0.0 },
  { 891.62834856284098, 363, 0.0 },
  { 894.80340366162693, 364, 0.0 },
  { 897.52824591642332, 365, 0.0 },
  { 900.70603699502817, 366, 0.0 },
  { 903.43360776447798, 367, 0.0 },
  { 906.61411993319712, 368, 0.0 },
  { 909.34440440851847, 369, 0.0 },
  { 912.52762293883552, 370, 0.0 },
  { 915.26060647112581, 371, 0.0 },
  { 918.44651679310869, 372, 0.0 },
  { 921.18218489076980, 373, 0.0 },
  { 924.37077259052307, 374, 0.0 },
  { 927.10911091674006, 375, 0.0 },
  { 930.30036173391295, 376, 0.0 },
  { 933.04135610418825, 377, 0.0 },
  { 936.23525592953251, 378, 0.0 },
  { 938.97889230927046, 379, 0.0 },
  { 942.17542718225297, 380, 0.0 },
  { 944.92169168439727, 381, 0.0 },
  { 948.12084779085956, 382, 0.0 },
  { 950.86972667357782, 383, 0.0 },
  { 954.07149034344741, 384, 0.0 },
  { 956.82297000786582, 385, 0.0 },
  { 960.02732771291210, 386, 0.0 },
  { 962.78139470089548, 387, 0.0 },
  { 965.98833305253549, 388, 0.0 },
  { 968.74497404451404, 389, 0.0 },
  { 971.95447979165920, 390, 0.0 },
  { 974.71368160449924, 391, 0.0 },
  { 977.92574163144968, 392, 0.0 },
  { 980.68749121636847, 393, 0.0 },
  { 983.90209254074750, 394, 0.0 },
  { 986.66637698126965, 395, 0.0 },
  { 989.88350675200218, 396, 0.0 },
  { 992.65031326195685, 397, 0.0 },
  { 995.86995875728655, 398, 0.0 },
  { 998.63927467884662, 399, 0.0 },
  { 1001.8614233043945, 400, 0.0 },
  { 1004.6332361061533, 401, 0.0 },
  { 1007.8578753930135, 402, 0.0 },
  { 1010.6321726681000, 403, 0.0 },
  { 1013.8592902709745, 404, 0.0 },
  { 1016.6360597352066, 405, 0.0 },
  { 1019.8656434305764, 406, 0.0 },
  { 1022.6448729206492, 407, 0.0 },
  { 1025.8769106049806, 408, 0.0 },
  { 1028.6585880766920, 409, 0.0 },
  { 1031.8930677646790, 410, 0.0 },
  { 1034.6771812911882, 411, 0.0 },
  { 1037.9140911140285, 412, 0.0 },
  { 1040.7006288841490, 413, 0.0 },
  { 1043.9399570878536, 414, 0.0 },
  { 1046.7289074043797, 415, 0.0 },
  { 1049.9706423481150, 416, 0.0 },
  { 1052.7619936261788, 417, 0.0 },
  { 1056.0061237806397, 418, 0.0 },
  { 1058.7998645461007, 419, 0.0 },
  { 1062.0463784919170, 420, 0.0 },
  { 1064.8424973797833, 421, 0.0 },
  { 1068.0913838059530, 422, 0.0 },
  { 1070.8898695588296, 423, 0.0 },
  { 1074.1411172611849, 424, 0.0 },
  { 1076.9419587277541, 425, 0.0 },
  { 1080.1955566074544, 426, 0.0 },
  { 1082.9987427409824, 427, 0.0 },
  { 1086.2546798030362, 428, 0.0 },
  { 1089.0601996599105, 429, 0.0 },
  { 1092.3184650117237, 430, 0.0 },
  { 1095.1263077500143, 431, 0.0 },
  { 1098.3868905999680, 432, 0.0 },
  { 1101.1970454780167, 433, 0.0 },
  { 1104.4599351340685, 434, 0.0 },
  { 1107.2723915091055, 435, 0.0 },
  { 1110.5375773774174, 436, 0.0 },
  { 1113.3523247042010, 437, 0.0 },
  { 1116.6197962877939, 438, 0.0 },
  { 1119.4368241172763, 439, 0.0 },
  { 1122.7065710147060, 440, 0.0 },
  { 1125.5258689927232, 441, 0.0 },
  { 1128.7978808967839, 442, 0.0 },
  { 1131.6194387627681, 443, 0.0 },
  { 1134.8937054592161, 444, 0.0 },
  { 1137.7175130449345, 445, 0.0 },
  { 1140.9940244112361, 446, 0.0 },
  { 1143.8200716395479, 447, 0.0 },
  { 1147.0988176436510, 448, 0.0 },
  { 1149.9270945272901, 449, 0.0 },
  { 1153.2080652264156, 450, 0.0 },
  { 1156.0385618667931, 451, 0.0 },
  { 1159.3217474062478, 452, 0.0 },
  { 1162.1544539922761, 453, 0.0 },
  { 1165.4398446042890, 454, 0.0 },
  { 1168.2747514112268, 455, 0.0 },
  { 1171.5623374138036, 456, 0.0 },
  { 1174.3994348021213, 457, 0.0 },
  { 1177.6892065979175, 458, 0.0 },
  { 1180.5284850121818, 459, 0.0 },
  { 1183.8204330874009, 460, 0.0 },
  { 1186.6618830551786, 461, 0.0 },
  { 1189.9559979784826, 462, 0.0 },
  { 1192.7996101092647, 463, 0.0 },
  { 1196.0958825307089, 464, 0.0 },
  { 1198.9416475148519, 465, 0.0 },
  { 1202.2400681648342, 466, 0.0 },
  { 1205.0879767725207, 467, 0.0 },
  { 1208.3885364607522, 468, 0.0 },
  { 1211.2385795409673, 469, 0.0 },
  { 1214.5412691554561, 470, 0.0 },
  { 1217.3934376349835, 471, 0.0 },
  { 1220.6982481410416, 472, 0.0 },
  { 1223.5525330234755, 473, 0.0 },
  { 1226.8594554627366, 474, 0.0 },
  { 1229.7158478275101, 475, 0.0 },
  { 1233.0248733169678, 476, 0.0 },
  { 1235.8833643183984, 477, 0.0 },
  { 1239.1944840494598, 478, 0.0 },
  { 1242.0550649158095, 479, 0.0 },
  { 1245.3682701533617, 480, 0.0 },
  { 1248.2309321859152, 481, 0.0 },
  { 1251.5462142674121, 482, 0.0 },
  { 1254.4109488395675, 483, 0.0 },
  { 1257.7282991741290, 484, 0.0 },
  { 1260.5950977305051, 485, 0.0 },
  { 1263.9145077980293, 486, 0.0 },
  { 1266.7833618535878, 487, 0.0 },
  { 1270.1048232038825, 488, 0.0 },
  { 1272.9757243430627, 489, 0.0 },
  { 1276.2992285949872, 490, 0.0 },
  { 1279.1721684708571, 491, 0.0 },
  { 1282.4977073114796, 492, 0.0 },
  { 1285.3726776448998, 493, 0.0 },
  { 1288.7002428286673, 494, 0.0 },
  { 1291.5772354074686, 495, 0.0 },
  { 1294.9068187553924, 496, 0.0 },
  { 1297.7858254335649, 497, 0.0 },
  { 1301.1174188324169, 498, 0.0 },
  { 1303.9984315293168, 499, 0.0 },
  { 1307.3320269308390, 500, 0.0 },
};
const double toler001 = 2.5000000000000017e-10;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_ldouble_factorial<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::ldouble_factorial<Ret>(data[i].n);
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
    VERIFY(!failure && max_abs_frac < toler);
    return num_errors;
  }

int
main()
{
  return test(data001, toler001);
}
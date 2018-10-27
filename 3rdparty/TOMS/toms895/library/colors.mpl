##### =========================
##### ===== COLORS MODULE =====
##### =========================

# in maple only 25 colors are predefined by default; this package
# provides 115 extra color names, as well as the predefined ones
# (for reference) - see http://www.learningwebdesign.com/colornames.html
# for samples of each color

colors := module()

  description "Colors package";

  export
    # predefined colors
    aquamarine, black, blue, brown, coral, cyan, gold, gray, green, khaki,
    magenta, maroon, navy, orange, pink, plum, red, sienna, tan, turquoise,
    violet, wheat, white, yellow,
    # new colors    
    aliceblue, antiquewhite, aqua, azure, beige, bisque, blanchedalmond,
    blueviolet, burlywood, cadetblue, chartreuse, chocolate, cornflowerblue,
    cornsilk, crimson, darkblue, darkcyan, darkgoldenrod, darkgray, darkgreen,
    darkkhaki, darkmagenta, darkolivegreen, darkorange, darkorchid, darkred,
    darksalmon, darkseagreen, darkslateblue, darkslategray, darkturquoise,
    darkviolet, deeppink, deepskyblue, dimgray, dodgerblue, firebrick,
    floralwhite, forestgreen, fuchsia, gainsboro, ghostwhite, goldenrod,
    greenyellow, honeydew, hotpink, indianred, indigo, ivory, lavender,
    lavenderblush, lawngreen, lemonchiffon, lighblue, lightcoral, lightcyan,
    lightgoldenrodyellow, lightgreen, lightgrey, lightpink, lightsalmon,
    lightseagreen, lightskyblue, lightslategray, lightsteelblue, lightyellow,
    lime, limegreen, linen, mediumaquamarine, mediumblue, mediumorchid,
    mediumpurple, mediumseagreen, mediumslateblue, mediumspringgreen,
    mediumturquoise, mediumvioletred, midnightblue, mintcream, mistyrose,
    moccasin, navajowhite, oldlace, olive, olivedrab, orangered, orchid,
    palegoldenrod, palegreen, paleturquoise, palevioletred, papayawhip,
    peachpuff, peru, powderblue, purple, rosybrown, royalblue, saddlebrown,
    salmon, sandybrown, seagreen, seashell, silver, skyblue, slateblue,
    slategray, snow, springgreen, steelblue, teal, thistle, tomato, whitesmoke,
    yellowgreen;
  
  option
    package;

# include the predefined colors for reference

  aquamarine := COLOR( RGB, 127/255, 255/255, 212/255 ):
  black := COLOR( RGB, 0/255, 0/255, 0/255 ):
  blue := COLOR( RGB, 0/255, 0/255, 255/255 ):
  brown := COLOR( RGB, 165/255, 42/255, 42/255 ):
  coral := COLOR( RGB, 255/255, 127/255, 80/255 ):
  cyan := COLOR( RGB, 0/255, 255/255, 255/255 ):
  gold := COLOR( RGB, 255/255, 215/255, 0/255 ):
  gray := COLOR( RGB, 128/255, 128/255, 128/255 ):
  green := COLOR( RGB, 0/255, 128/255, 0/255 ):
  khaki := COLOR( RGB, 240/255, 230/255, 140/255 ):
  magenta := COLOR( RGB, 255/255, 0/255, 255/255 ):
  maroon := COLOR( RGB, 128/255, 0/255, 0/255 ):
  navy := COLOR( RGB, 0/255, 0/255, 128/255 ):
  orange := COLOR( RGB, 255/255, 165/255, 0/255 ):
  pink := COLOR( RGB, 255/255, 192/255, 203/255 ):
  plum := COLOR( RGB, 221/255, 160/255, 221/255 ):
  red := COLOR( RGB, 225/255, 0/255, 0/255 ):
  sienna := COLOR( RGB, 160/255, 82/255, 45/255 ):
  tan := COLOR( RGB, 210/255, 180/255, 140/255 ):
  turquoise := COLOR( RGB, 64/255, 224/255, 208/255 ):
  violet := COLOR( RGB, 238/255, 130/255, 238/255 ):
  wheat := COLOR( RGB, 245/255, 222/255, 179/255 ):
  white := COLOR( RGB, 255/255, 255/255, 255/255 ):
  yellow := COLOR( RGB, 255/255, 255/255, 0/255 ):
    
# NEW COLOR DEFINITIONS

  aliceblue := COLOR( RGB, 240/255, 248/255, 255/255 ):
  antiquewhite := COLOR( RGB, 250/255, 235/255, 215/255 ):
  aqua := COLOR( RGB, 0/255, 255/255, 255/255 ):
  azure := COLOR( RGB, 240/255, 255/255, 255/255 ):
  beige := COLOR( RGB, 245/255, 245/255, 220/255 ):
  bisque := COLOR( RGB, 255/255, 228/255, 196/255 ):
  blanchedalmond := COLOR( RGB, 255/255, 255/255, 205/255 ):
  blueviolet := COLOR( RGB, 138/255, 43/255, 226/255 ):
  burlywood := COLOR( RGB, 222/255, 184/255, 135/255 ):
  cadetblue := COLOR( RGB, 95/255, 158/255, 160/255 ):
  chartreuse := COLOR( RGB, 127/255, 255/255, 0/255 ):
  chocolate := COLOR( RGB, 210/255, 105/255, 30/255 ):
  cornflowerblue := COLOR( RGB, 100/255, 149/255, 237/255 ):
  cornsilk := COLOR( RGB, 255/255, 248/255, 220/255 ):
  crimson := COLOR( RGB, 220/255, 20/255, 60/255 ):
  darkblue := COLOR( RGB, 0/255, 0/255, 139/255 ):
  darkcyan := COLOR( RGB, 0/255, 139/255, 139/255 ):
  darkgoldenrod := COLOR( RGB, 184/255, 134/255, 11/255 ):
  darkgray := COLOR( RGB, 169/255, 169/255, 169/255 ):
  darkgreen := COLOR( RGB, 0/255, 100/255, 0/255 ):
  darkkhaki := COLOR( RGB, 189/255, 183/255, 107/255 ):
  darkmagenta := COLOR( RGB, 139/255, 0/255, 139/255 ):
  darkolivegreen := COLOR( RGB, 85/255, 107/255, 47/255 ):
  darkorange := COLOR( RGB, 255/255, 140/255, 0/255 ):
  darkorchid := COLOR( RGB, 153/255, 50/255, 204/255 ):
  darkred := COLOR( RGB, 139/255, 0/255, 0/255 ):
  darksalmon := COLOR( RGB, 233/255, 150/255, 122/255 ):
  darkseagreen := COLOR( RGB, 143/255, 188/255, 143/255 ):
  darkslateblue := COLOR( RGB, 72/255, 61/255, 139/255 ):
  darkslategray := COLOR( RGB, 47/255, 79/255, 79/255 ):
  darkturquoise := COLOR( RGB, 0/255, 206/255, 209/255 ):
  darkviolet := COLOR( RGB, 148/255, 0/255, 211/255 ):
  deeppink := COLOR( RGB, 255/255, 20/255, 147/255 ):
  deepskyblue := COLOR( RGB, 0/255, 191/255, 255/255 ):
  dimgray := COLOR( RGB, 105/255, 105/255, 105/255 ):
  dodgerblue := COLOR( RGB, 30/255, 144/255, 255/255 ):
  firebrick := COLOR( RGB, 178/255, 34/255, 34/255 ):
  floralwhite := COLOR( RGB, 255/255, 250/255, 240/255 ):
  forestgreen := COLOR( RGB, 34/255, 139/255, 34/255 ):
  fuchsia := COLOR( RGB, 255/255, 0/255, 255/255 ):
  gainsboro := COLOR( RGB, 220/255, 220/255, 220/255 ):
  ghostwhite := COLOR( RGB, 248/255, 248/255, 255/255 ):
  goldenrod := COLOR( RGB, 218/255, 165/255, 32/255 ):
  greenyellow := COLOR( RGB, 173/255, 255/255, 47/255 ):
  honeydew := COLOR( RGB, 240/255, 255/255, 240/255 ):
  hotpink := COLOR( RGB, 255/255, 105/255, 180/255 ):
  indianred := COLOR( RGB, 205/255, 92/255, 92/255 ):
  indigo := COLOR( RGB, 75/255, 0/255, 130/255 ):
  ivory := COLOR( RGB, 255/255, 240/255, 240/255 ):
  lavender := COLOR( RGB, 230/255, 230/255, 250/255 ):
  lavenderblush := COLOR( RGB, 255/255, 240/255, 245/255 ):
  lawngreen := COLOR( RGB, 124/255, 252/255, 0/255 ):
  lemonchiffon := COLOR( RGB, 255/255, 250/255, 205/255 ):
  lighblue := COLOR( RGB, 173/255, 216/255, 230/255 ):
  lightcoral := COLOR( RGB, 240/255, 128/255, 128/255 ):
  lightcyan := COLOR( RGB, 224/255, 255/255, 255/255 ):
  lightgoldenrodyellow := COLOR( RGB, 250/255, 250/255, 210/255 ):
  lightgreen := COLOR( RGB, 144/255, 238/255, 144/255 ):
  lightgrey := COLOR( RGB, 211/255, 211/255, 211/255 ):
  lightpink := COLOR( RGB, 255/255, 182/255, 193/255 ):
  lightsalmon := COLOR( RGB, 255/255, 160/255, 122/255 ):
  lightseagreen := COLOR( RGB, 32/255, 178/255, 170/255 ):
  lightskyblue := COLOR( RGB, 135/255, 206/255, 250/255 ):
  lightslategray := COLOR( RGB, 119/255, 136/255, 153/255 ):
  lightsteelblue := COLOR( RGB, 176/255, 196/255, 222/255 ):
  lightyellow := COLOR( RGB, 255/255, 255/255, 224/255 ):
  lime := COLOR( RGB, 0/255, 255/255, 0/255 ):
  limegreen := COLOR( RGB, 50/255, 205/255, 50/255 ):
  linen := COLOR( RGB, 250/255, 240/255, 230/255 ):
  mediumaquamarine := COLOR( RGB, 102/255, 205/255, 170/255 ):
  mediumblue := COLOR( RGB, 0/255, 0/255, 205/255 ):
  mediumorchid := COLOR( RGB, 186/255, 85/255, 211/255 ):
  mediumpurple := COLOR( RGB, 147/255, 112/255, 219/255 ):
  mediumseagreen := COLOR( RGB, 60/255, 179/255, 113/255 ):
  mediumslateblue := COLOR( RGB, 123/255, 104/255, 238/255 ):
  mediumspringgreen := COLOR( RGB, 0/255, 250/255, 154/255 ):
  mediumturquoise := COLOR( RGB, 72/255, 209/255, 204/255 ):
  mediumvioletred := COLOR( RGB, 199/255, 21/255, 133/255 ):
  midnightblue := COLOR( RGB, 25/255, 25/255, 112/255 ):
  mintcream := COLOR( RGB, 245/255, 255/255, 250/255 ):
  mistyrose := COLOR( RGB, 255/255, 228/255, 225/255 ):
  moccasin := COLOR( RGB, 255/255, 228/255, 181/255 ):
  navajowhite := COLOR( RGB, 255/255, 222/255, 173/255 ):
  oldlace := COLOR( RGB, 253/255, 245/255, 230/255 ):
  olive := COLOR( RGB, 128/255, 128/255, 0/255 ):
  olivedrab := COLOR( RGB, 107/255, 142/255, 35/255 ):
  orangered := COLOR( RGB, 255/255, 69/255, 0/255 ):
  orchid := COLOR( RGB, 218/255, 112/255, 214/255 ):
  palegoldenrod := COLOR( RGB, 238/255, 232/255, 170/255 ):
  palegreen := COLOR( RGB, 152/255, 251/255, 152/255 ):
  paleturquoise := COLOR( RGB, 175/255, 238/255, 238/255 ):
  palevioletred := COLOR( RGB, 219/255, 112/255, 147/255 ):
  papayawhip := COLOR( RGB, 255/255, 239/255, 213/255 ):
  peachpuff := COLOR( RGB, 255/255, 239/255, 213/255 ):
  peru := COLOR( RGB, 205/255, 133/255, 63/255 ):
  powderblue := COLOR( RGB, 176/255, 224/255, 230/255 ):
  purple := COLOR( RGB, 128/255, 0/255, 128/255 ):
  rosybrown := COLOR( RGB, 188/255, 143/255, 143/255 ):
  royalblue := COLOR( RGB, 65/255, 105/255, 225/255 ):
  saddlebrown := COLOR( RGB, 139/255, 69/255, 19/255 ):
  salmon := COLOR( RGB, 250/255, 128/255, 114/255 ):
  sandybrown := COLOR( RGB, 244/255, 164/255, 96/255 ):
  seagreen := COLOR( RGB, 46/255, 139/255, 87/255 ):
  seashell := COLOR( RGB, 255/255, 245/255, 238/255 ):
  silver := COLOR( RGB, 192/255, 192/255, 192/255 ):
  skyblue := COLOR( RGB, 135/255, 206/255, 235/255 ):
  slateblue := COLOR( RGB, 106/255, 90/255, 205/255 ):
  slategray := COLOR( RGB, 112/255, 128/255, 144/255 ):
  snow := COLOR( RGB, 255/255, 250/255, 250/255 ):
  springgreen := COLOR( RGB, 0/255, 255/255, 127/255 ):
  steelblue := COLOR( RGB, 70/255, 130/255, 180/255 ):
  teal := COLOR( RGB, 0/255, 128/255, 128/255 ):
  thistle := COLOR( RGB, 216/255, 191/255, 216/255 ):
  tomato := COLOR( RGB, 253/255, 99/255, 71/255 ):
  whitesmoke := COLOR( RGB, 245/255, 245/255, 245/255 ):
  yellowgreen := COLOR( RGB, 154/255, 205/255, 50/255 ):

end module:

##### ===== END OF COLORS MODULE =====

# add the colors module to the CFSF.lib repository

try
  savelibname := cat( currentdir(), "/CFSF.lib" ):
  try march( 'create', savelibname, 500 ) catch: end try:
  savelib( 'colors' ):
  map( march, ['gc','reindex','pack'], savelibname ):
end try:

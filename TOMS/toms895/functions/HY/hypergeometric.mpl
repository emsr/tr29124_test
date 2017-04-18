create( 'series',
  label = "HY.hyper.power.01",
  booklabel = "15.1.4",
  general = [ pochhammer(a, k) * pochhammer(b, k) * z^k / (pochhammer(c, k) * k!) ],
  parameters = { a, b, c },
  constraints = { c::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a, b], [c], z),
  category = "power series"
):

create( 'contfrac',
  label = "HY.hyper.cfrac.01",
  booklabel = "15.3.3",
  front = 1,
  even = [-( ( b + ( m / 2 ) ) * ( c - a + ( m / 2 ) ) ) / ( ( c + m - 1 ) * ( c + m ) ) * z, 1],
  odd = [-( a + ( ( m - 1 ) / 2 ) ) * ( c - b + ( ( m - 1 ) / 2 ) ) / ( ( c + m - 1 ) * ( c + m ) ) * z, 1],
  parameters = {a, b, c},
  constraints = { z::Not(RealRange(1,infinity)), c::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a, b], [c], z) / hypergeom([a, b + 1], [c + 1], z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "HY.hyper.cfrac.02",
  booklabel = "15.3.4",
  begin = [[z, 1]],
  even = [-( ( ( a + ( m - 2 )/ 2 ) * ( c + ( ( m - 2 ) / 2 ) ) ) / ( ( c + m - 2 ) * ( c +  m - 1 ) ) ) * z , 1],
  odd = [-( ( ( ( m - 1 ) / 2 ) * ( c - a + ( ( m - 1 ) / 2 ) ) ) / ( ( c + m - 2 ) * ( c + m - 1 ) ) ) * z, 1],
  parameters = {a, c},
  constraints = { z::Not(RealRange(1,infinity)), c::Not(nonposint) },
  function = hypergeom,
  lhs = z * hypergeom([a, 1], [c + 1],  z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "HY.hyper.cfrac.03",
  booklabel = "15.3.7",
  factor = 1/z,
  begin = [[z, 1]],
  general = [[-(m-1)^2 * z / (4*(m-1)^2 - 1), 1]],
  constraints = { z::Not(RealRange(1,infinity)) },
  function = hypergeom,
  lhs = hypergeom([1/2, 1], [3/2], z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "HY.hyper.tfrac.01",
  booklabel = "15.3.8",
  front = (c + ( b - a + 1 ) * z) / c,
  general = [ [-( c - a + m ) * ( b + m ) * z, c + m + ( b - a + m + 1 ) * z] ],
  parameters = {a, b, c},
  constraints = { abs(z) < 1, c::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a, b], [c], z) / hypergeom([a, b + 1], [c + 1], z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "HY.hyper.mfrac.01",
  booklabel = "15.3.9a",
  begin = [[ c, c+(1-a)*z ]],
  general = [ [-(m-1) * z * (c-a+m-1), c+m-1+(m-a)*z ] ],
  parameters = {a, c},
  constraints = { abs(z) < 1, c::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a, 1], [c+1], z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "HY.hyper.mfrac.02",
  booklabel = "15.3.9b",
  factor = (1-a) * z / c,
  begin = [[ c, c+(1-a)*z ]],
  general = [ [-(m-1) * z * (c-a+m-1), c+m-1+(m-a)*z ] ],
  parameters = {a, c},
  constraints = { abs(z) < 1, a::Not(nonposint), a <> 1, a <> 2 },
  function = hypergeom,
  lhs = hypergeom([1-c, 1], [2-a], 1/z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "HY.hyper.mfrac.03",
  booklabel = "15.3.12",
  begin = [[ 1/2, 1/2+1/2*z ]],
  general = [ [-(m-1) * z * (m-1), 1/2+m-1+(m-1/2)*z ] ],
  constraints = { abs(z) < 1 },
  function = hypergeom,
  lhs = hypergeom([1/2, 1], [3/2], z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "HY.hyper.tfrac.02",
  booklabel = "15.3.13",
  front = 1 - ( ( a + b + 1 ) / c ) * z,
  general = [ [( ( ( a + m ) * ( b + m ) ) / ( ( c + m - 1 ) * ( c + m ) ) ) * ( z - z^2 ), 1 - ( ( a + b + 2*m + 1 ) / ( c + m ) ) * z] ],
  parameters = {a, b, c},
  constraints = { Re(z) < 1/2, c::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a, b], [c], z) / hypergeom([a + 1, b + 1], [c + 1], z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "HY.hyper.mfrac.04",
  booklabel = "15.3.14",
  begin = [[c, c-(a+1)*z]],
  general = [[(a+m)*m * (z-z^2), c+m - (a+2*m+1)*z ]],
  parameters = {a, c},
  function = hypergeom,
  lhs = hypergeom([a+1, 1], [c+1], z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "HY.hyper.norlund.01",
  booklabel = "15.3.17",
  begin = [[1, 1 - z]],
  general = [ [ (m-1)/(m-1 +1/2) * (z - z^2), 1 - ( 2*(m-1) +1/2 )/(m-1 +1/2) * z ] ],
  constraints = { Re(z) < 1/2 },
  function = hypergeom,
  lhs = hypergeom([1/2, 1], [3/2], z),
  category = "Norlund fraction"
):

create( 'contfrac',
  label = "HY.hyper.frac.01",
  booklabel = "15.6.4",
  front = 1,
  begin = [[ -b*c/d, e-a-1 ]],
  odd = [ - (d - a + (m-1)/2 - 1 ) * (b + (m-1)/2) * (c + (m-1)/2) / ( (d + m-2 ) * (d + m-1) ), e - a - 1 ],
  even = [ (a + m/2) * (d - b + m/2 - 1) * (d - c + m/2 - 1) / ( (d + m-2 ) * (d + m-1) ), 1 ],
  parameters = {a, b, c, d, e},
  function = hypergeom,
  lhs = hypergeom([a, b, c], [d, e], 1) / hypergeom([a+1, b, c], [d, e], 1),
  category = "continued fraction"
):

create( 'contfrac',
  label = "HY.hyper.frac.02",
  booklabel = "15.6.5",
  front = 1,
  odd = [( ( a + ( m - 1 ) / 2 ) * ( b + ( m - 1 ) / 2 ) * ( c + ( m - 1 ) / 2 ) ) / ( ( d + m - 1  ) * ( d + m  ) ),  d + e - a - b - c],
  even = [ -( ( ( m / 2 ) + d - a ) * ( ( m / 2 ) + d - b ) * ( ( m / 2 ) + d - c ) ) / ( ( d + m - 1 ) * ( d + m  ) ), 1],
  parameters = {a, b, c, d, e},
  function = hypergeom,
  lhs = hypergeom([a, b, c], [d, e], 1) / hypergeom([a, b, c], [d + 1, e], 1),
  category = "continued fraction"
):

create( 'contfrac',
  label = "HY.hyper.frac.03",
  booklabel = "15.6.6",
  front = 1,
  odd = [-( ( ( m - 1 ) / 2 + d - a ) * ( ( m - 1 ) / 2 + b ) * ( ( m - 1 ) / 2 + c ) ) / ( ( d + m - 1 ) * ( d + m) ), (e - a - 1)],
  even = [( ( a + m / 2 ) * ( d - b + m / 2 ) * ( d - c + m / 2 ) ) / ( ( d + m  - 1 ) * ( d + m  ) ), 1],
  parameters = {a, b, c, d, e},
  function = hypergeom,
  lhs = hypergeom([a, b, c], [d, e], 1) / hypergeom([a + 1, b, c], [d + 1, e], 1),
  category = "continued fraction"
):

create( 'contfrac',
  label = "HY.hyper.frac.04",
  booklabel = "15.6.7",
  front = (e - (a-1) - 1) / e,
  odd = [((a-1)+(m+1)/2) * ((d-1)-b+(m+1)/2) * ((d-1)-c+(m+1)/2) / ( (e+(m-1)/2) * ((d-1)+m) * ((d-1)+m+1) ) , 1],
  even = [- ((d-1)-(a-1)+m/2) * (b+m/2) * (c+m/2) / ( (e+(m/2)) * ((d-1)+m) *((d-1)+m+1) ) , (e-(a-1)-1)/(e+m/2)],
  parameters = {a, b, c, d, e},
  function = hypergeom,
  lhs = hypergeom([a, b, c], [d, e], 1) / hypergeom([a, b+1, c+1], [d+1, e+1], 1),
  category = "continued fraction"
):

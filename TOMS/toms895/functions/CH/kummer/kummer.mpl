create( 'series',
  label = "CH.confl.power.01",
  booklabel = "16.1.2",
  general = [pochhammer(a, k) / pochhammer(b, k) * z^k / k! ],
  parameters = {a, b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a], [b], z),
  category = "power series"
);

create( 'series',
  label = "CH.confl.power.02",
  booklabel = "16.1.12",
  general = [pochhammer(a, k) * pochhammer(b, k) * z^k / k! ],
  parameters = {a, b},
  function = hypergeom,
  lhs = hypergeom([a, b], [], z),
  category = "power series"
);

create( 'contfrac',
  label = "CH.confl.cfrac.01",
  booklabel = "16.1.13",
  front = 1,
  even = [z* ( a + m / 2 ) / ( ( b + m - 1 ) * ( b + m ) ), 1],
  odd = [-z * ( b - a + ( ( m - 1 ) / 2 ) ) / ( ( b + m - 1 ) * ( b + m ) ), 1],
  parameters = {a, b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a], [b], z) / hypergeom([a + 1], [b + 1], z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "CH.confl.cfrac.02",
  booklabel = "16.1.14",
  begin = [[z, 1]],
  even = [-z * ( b + m/2 -1 ) / ( ( b + m - 2 ) * ( b + m - 1 ) ), 1],
  odd = [z * (  ( m - 1 ) / 2 ) / ( ( b + m - 2 ) * ( b + m - 1 ) ), 1],
  parameters = {b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = z * hypergeom([1], [b+1], z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "CH.confl.tfrac.01",
  booklabel = "16.1.16",
  front = (b-z) / b,
  factor = 1/b,
  general = [[(a+m)*z, b+m-z ]],
  parameters = {a, b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([a], [b], z) / hypergeom([a+1], [b+1], z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "CH.confl.mfrac.01",
  booklabel = "16.1.17",
  begin = [[b, b - z]],
  general = [[(m-1)*z, b+m-1-z ]],
  parameters = {b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([1], [b+1], z),
  category = "M-fraction"
):


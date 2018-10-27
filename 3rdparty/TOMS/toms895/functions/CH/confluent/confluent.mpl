create( 'contfrac',
  label = "CH.confl.cfrac.03",
  booklabel = "16.2.4",
  front = 1,
  even = [-( b + m / 2 ) * z, 1],
  odd = [-( a + ( ( m - 1 ) / 2) ) * z, 1],
  parameters = {a, b},
  constraints = { z::Not(positive) },
  function = hypergeom,
  lhs = hypergeom([a, b], [], z) / hypergeom([a, b + 1], [], z),
  category = "C-fraction"
):


create( 'series',
  label = "CH.confl.power.03",
  booklabel = "16.3.1",
  general = [ z^k / (pochhammer(b, k) * k!) ],
  parameters = {b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([], [b], z),
  category = "power series"

):

create( 'contfrac',
  label = "CH.confl.cfrac.04",
  booklabel = "16.3.4",
  front = 1,
  general = [ [1 / ( ( b + m - 1 ) * ( b + m ) ) * z, 1] ],
  parameters = {b},
  constraints = { b::Not(nonposint) },
  function = hypergeom,
  lhs = hypergeom([], [b], z) / hypergeom([], [b + 1], z),
  category = "C-fraction"
):

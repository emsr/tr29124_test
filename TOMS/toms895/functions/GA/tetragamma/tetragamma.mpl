create( 'contfrac',
  label = "GA.psi2.sfrac.01",
  booklabel = "12.5.1",
  front = -1 / z^2 - 1 / z^3,
  factor = -( ( 2 * Pi ) / z)^2,
  begin = [[( 1 / ( 8 * Pi^2 ) ), z^2]],
  even = [( ( m / 2 )^2 * ( ( m / 2 ) + 1 ) ) / (2 * ( m + 1 ) ), 1],
  odd = [( ( ( m - 1 ) / 2 ) * ( ( ( m - 1 ) / 2 ) + 1 )^2 ) / (2 * m ), z^2],
  constraints = { abs(functions:-argument(z)) < Pi/2 },
  function = Psi,
  lhs = Psi(2, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "GA.psi2.sfrac.02",
  booklabel = "12.5.2",
  factor = -1,
  begin = [[1 / ( z * ( z - 1 ) ), 1]],
  even = [( m / 2 )^4 / ( z * (z - 1) ), m],
  odd = [( ( m - 1)/2)^4 / ( z * (z - 1) ), m],
  constraints = { Re(z) > 1/2, z::Not(RealRange(Open(1/2),1)) },
  function = Psi,
  lhs = Psi(2, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "GA.psi2.cfrac.01",
  booklabel = "12.5.3",
  front = -1/z,
  begin = [[1/z,1]],
  general = [[( ( ( m + 2 ) / 4 )^2 - 2 * ( ( m + 2 ) / 4 ) + 2 ) / ( 2 * ( ( m + 2 ) / 4 ) - 1 ) * z^(-1), 1], [-( ( ( m + 1 ) / 4 )^2 + 1 ) / ( 2 * ( ( m + 1 ) / 4 ) - 1 ) * z^(-1), 1], [( m / 4 )^3 / ( 2 * ( ( m / 4 )^2 + 1 ) ) * z^(-1), 1], [-( ( m - 1 ) / 4 )^3 / ( 2 * ( ( ( m - 1 ) / 4 )^2 + 1 ) ) * z^(-1), 1]],
  constraints = { Re(z) > 1 },
  function = Psi,
  lhs = Psi(2, z),
  category = "C-fraction"
):

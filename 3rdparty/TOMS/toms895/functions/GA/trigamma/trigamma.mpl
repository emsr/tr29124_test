create( 'contfrac',
  label = "GA.psi1.sfrac.01",
  booklabel = "12.4.1",
  front = 1 / z + 1 / ( 2 * z^2 ),
  factor = 2 * Pi / z,
  begin = [[(1 / (12 * Pi)), z^2]],
  even = [( m^2 * ( m^2 - 1 ) ) / (4 * ( 4 * m^2 -1 ) ), 1],
  odd = [( m^2 * ( m^2 - 1 ) ) / (4 * ( 4 * m^2 -1 ) ), z^2],
  constraints = { abs(functions:-argument(z)) < Pi/2 },
  function = Psi,
  lhs = Psi(1, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "GA.psi1.cfrac.01",
  booklabel = "12.4.2",
  begin = [[z^( -1 ) , 1]],
  even = [( -( m / 2)^2 / ( 2 * m - 2) ) * z^( -1 ), 1],
  odd = [( ( ( m - 1 ) / 2 )^2 / ( 2 * m) ) * z^( -1 ), 1],
  constraints = { Re(z) > 1/2 },
  function = Psi,
  lhs = Psi(1, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "GA.psi1.jfrac.01",
  booklabel = "12.4.3",
  begin = [[1 , -1/2 + z]],
  general = [ [( (m - 1 )^4 ) / ( 4 * ( 2 * m-3 ) * ( 2 * m-1 ) ), -1/2 + z] ],
  constraints = { Re(z) > 1/2 },
  function = Psi,
  lhs = Psi(1, z),
  category = "J-fraction"
):

create( 'series',
  label = "GA.psi0.power.01",
  booklabel = "12.3.2a",
  front = -gamma,
  general = [ ( 1 / ( k ) ) - ( 1 / ( z + k - 1) ) ],
  constraints = { z::Not(nonposint) },
  function = Psi,
  lhs = Psi(z),
  category = "power series"
):

create( 'series',
  label = "GA.psi.power.01",
  booklabel = "12.3.2b",
  factor = ( -1 )^( alpha + 1 ) * alpha!,
  general = [ 1 / ( ( z + k )^( alpha + 1 ) ) ],
  parameters = { alpha },
  constraints = { z::Not(nonposint), alpha::posint },
  function = Psi,
  lhs = Psi(alpha, z),
  category = "power series"
):

create( 'series',
  label = "GA.psi0.asymp.01",
  booklabel = "12.3.7",
  front = ln(z) - ( 1 / ( 2 * z ) ),
  general = [ -( bernoulli( 2 * k ) ) / ( 2 * k ) * z^( -2 * k ) ],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Psi,
  lhs = Psi(z),
  category = "asymptotic series"
):

create( 'series',
  label = "GA.psi.asymp.01",
  booklabel = "12.3.8",
  factor = ( -1 )^( alpha - 1 ),
  front = ( ( ( alpha - 1 )! / z^alpha ) + ( alpha! / ( 2 * z^( alpha + 1 ) ) ) ) * ( -1 )^( alpha - 1 ),
  general = [ ( ( bernoulli( 2 * k ) * ( 2 * k + alpha - 1 )! )  / ( 2 * k )! ) * z^( -2 * k - alpha) ],
  parameters = { alpha },
  constraints = { abs(functions:-argument(z)) < Pi, alpha::posint },
  function = Psi,
  lhs = Psi(alpha, z),
  category = "asymptotic series"
):

$ifdef QD
create( 'contfrac',
  label = "GA.psi0.sfrac.01",
  booklabel = "12.3.11",
  front = ln(z) - 1/(2*z),
  factor = (-1),
  even = [qdpolygamma(m, 0), 1],
  odd = [qdpolygamma(m, 0), z^2],
  parameters = { alpha },
  function = Psi,
  lhs = Psi(alpha, z),
  category = "S-fraction"
):
$endif

$ifdef QD
create( 'contfrac',
  label = "GA.psi.sfrac.01",
  booklabel = "12.3.12",
  front = (-1)^(alpha-1) * ((alpha-1)!/z^alpha + alpha!/(2*z^(alpha+1))),
  factor = (-1)^(alpha-1) * (2*Pi / z)^alpha,
  even = [qdpolygamma(m, alpha), 1],
  odd = [qdpolygamma(m, alpha), z^2],
  parameters = { alpha },
  constraints = { alpha::posint },
  function = Psi,
  lhs = Psi(alpha, z),
  category = "S-fraction"
):
$endif

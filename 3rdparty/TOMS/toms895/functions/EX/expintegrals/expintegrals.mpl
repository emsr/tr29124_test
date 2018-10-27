create( 'series',
  label = "EX.ei.power.01",
  booklabel = "14.1.10",
  front = -gamma - ln( z ),
  factor = -1,
  general = [ (-1)^k * z^k / ( k * k! ) ],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Ei,
  lhs = Ei(z),
  category = "power series"
):

create( 'series',
  label = "EX.ei.power.02",
  booklabel = "14.1.11",
  front = ( ( -z )^( n - 1 ) / ( n - 1 )! ) * ( -gamma -ln(z) + sum(1 / i, i = 1..n-1 ) ) - sum( ( -z )^i / ( (i - n + 1) * i! ) ,i=0..n-2),
  general = [ -( ( -z )^( k + n - 1) / ( ( k ) * ( k + n - 1)! ) ) ],
  parameters = { n },
  constraints = { abs(functions:-argument(z)) < Pi, n::posint },
  function = Ei,
  lhs = Ei(n, z),
  category = "power series"
):

create( 'series',
  label = "EX.ei.power.03",
  booklabel = "14.1.12",
  front = GAMMA(1 - nu) * z^(nu - 1) - 1/(1-nu),
  general = [ (-1)^(k+1) * z^k  / ( (k!) * ( k+1 -nu )) ],
  parameters = { nu },
  constraints = { z <> 0, nu::Not(posint) },
  function = Ei,
  lhs = Ei(nu, z),
  category = "power series"
):

create( 'series',
  label = "EX.ei.asymp.01",
  booklabel = "14.1.13",
  factor = exp(-z),
  general = [ ( (-1)^k * pochhammer(nu, k) ) / z^( k + 1 ) ],
  parameters = { nu },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Ei,
  lhs = Ei(nu, z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "EX.ei.sfrac.01",
  booklabel = "14.1.16",
  factor = exp(-z),
  begin = [[1, z],[n, 1]],
  even = [n + ( m - 2 ) / 2, 1],
  odd = [( m - 1 ) / 2, z],
  parameters = { n },
  constraints = { abs(functions:-argument(z)) < Pi, n::posint },
  function = Ei,
  lhs = Ei(n, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EX.ei.cfrac.01",
  booklabel = "14.1.19",
  begin = [[ exp(-z) / z, 1 ]],
  even = [ (m/2 + nu - 1) / z, 1 ],
  odd = [ ((m-1)/2) / z, 1 ],
  parameters = { nu },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Ei,
  lhs = Ei(nu, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "EX.ei.cfrac.02",
  booklabel = "14.1.20",
  front = z^( nu - 1) * GAMMA(1 - nu),
  begin = [[-exp(-z), 1 - nu]],
  even = [( nu - m / 2 ) * z, m - nu],
  odd = [( ( m - 1 ) / 2 ) * z, m - nu],
  parameters = { nu },
  constraints = { z <> 0, nu::Not(posint) },
  function = Ei,
  lhs = Ei(nu, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "EX.ei.mfrac.01",
  booklabel = "14.1.22",
  front = z^( nu - 1) * GAMMA(1 - nu),
  factor = -exp(-z),
  begin = [[1, 1 - nu - z]],
  general = [ [( m - 1 ) * z, m - nu - z] ],
  parameters = { nu },
  constraints = { z <> 0, nu::Not(posint) },
  function = Ei,
  lhs = Ei(nu, z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "EX.ei.jfrac.01",
  booklabel = "14.1.23",
  factor = exp(-z),
  begin = [[1, nu + z]],
  general = [ [-( m - 1 ) * ( nu + m - 2), nu + ( m - 1 )*2 + z] ],
  parameters = { nu },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Ei,
  lhs = Ei(nu, z),
  category = "J-fraction"
):

create( 'contfrac',
  label = "EX.ei.jfrac.02",
  booklabel = "14.1.24",
  front = exp(-z)/z,
  begin = [[-exp(-z) * nu, z * (1 + z) + nu * z]],
  general = [[-(m - 1) * (nu + m - 1) * z^2 , z * (m + z) + (nu + m - 1) * z]],
  parameters = { nu },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = Ei,
  lhs = Ei(nu, z),
  category = "J-fraction"
):

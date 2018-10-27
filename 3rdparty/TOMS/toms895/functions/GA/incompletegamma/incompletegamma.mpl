create( 'series',
  label = "GA.lincgamma.power.01",
  booklabel = "12.6.7",
  factor = z^a,
  general = [ ( -z )^k / ( ( a + k ) * k! ) ],
  parameters = { a },
  constraints = { Re(a) > 0 },
  function = { functions:-lincgamma, GAMMA },
  lhs = functions:-lincgamma(a, z),
  category = "power series"
):

create( 'series',
  label = "GA.lincgamma.power.02",
  booklabel = "12.6.8",
  factor = ( z^a * exp(-z) ) / a,
  general = [ z^k / pochhammer(1 + a, k) ],
  parameters = { a },
  constraints = { Re(a) > 0 },
  function = { functions:-lincgamma, GAMMA },
  lhs = functions:-lincgamma(a, z),
  category = "power series"
):

create( 'series',
  label = "GA.uincgamma.asymp.01",
  booklabel = "12.6.10",
  factor = z^a * exp(-z),
  general = [ ( -1 )^k * pochhammer(1 - a, k) * z^(-k-1) ],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "GA.uincgamma.sfrac.01",
  booklabel = "12.6.15",
  factor = z^a * exp(-z),
  begin = [[1, z]],
  even = [ ( m / 2 ) - a, 1 ],
  odd = [ ( m - 1 ) / 2 , z ],
  parameters = { a },
  constraints = { a::RealRange(-infinity, Open(1)), abs(functions:-argument(z)) < Pi },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "GA.uincgamma.cfrac.01",
  booklabel = "12.6.17",
  factor = z^a * exp(-z),
  begin = [[1/z, 1]],
  even = [ (( m / 2 ) - a) * z^(-1), 1 ],
  odd = [ (( m - 1 ) / 2) * z^(-1), 1 ],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "GA.lincgamma.cfrac.01",
  booklabel = "12.6.23",
  factor = z^( a-1 ) * exp( -z ),
  begin = [[z/a, 1]],
  even = [( -( a + ( m / 2 ) - 1 ) / ( ( a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
  odd = [ ( ( m - 1 ) / ( 2 * (a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
  parameters = { a },
  constraints = { Re(a) > 0 },
  function = { functions:-lincgamma, GAMMA },
  lhs = functions:-lincgamma(a, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "GA.uincgamma.cfrac.02",
  booklabel = "12.6.24",
  front = GAMMA( a ),
  factor = -( z^( a ) * exp( -z ) ) / ( z ),
  begin = [[z/a, 1]],
  even = [( -( a + ( m / 2 ) - 1 ) / ( ( a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
  odd = [ ( ( m - 1 ) / ( 2 * (a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
#  begin = [[z, 1], [( -1 / ( a + 1) ) * z, 1]],
#  even = [( -( a + ( m / 2 ) - 1 ) / ( ( a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
#  odd = [ ( ( m - 1 ) / ( 2 * (a + m - 2 ) * ( a + m - 1 ) ) ) * z, 1],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi, a::Not(nonposint) },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "GA.lincgamma.mfrac.01",
  booklabel = "12.6.30",
  factor = z^( a ) * exp( -z ),
  begin = [[ 1, a - z ]],
  general = [[ ( m - 1 ) * z, a + ( m - 1 ) - z ]],
  parameters = { a },
  constraints = { a::Not(nonposint) },
  function = { functions:-lincgamma, GAMMA },
  lhs = functions:-lincgamma(a, z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "GA.uincgamma.mfrac.01",
  booklabel = "12.6.31",
  front = GAMMA( a ),
  factor = z^( a ) * exp( -z ),
  begin = [[ -1, a - z ]],
  general = [[ (m - 1) * z, a + (m - 1) - z ]],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi, Re(a) > 0 },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "M-fraction"
):

create( 'contfrac',
  label = "GA.uincgamma.jfrac.01",
  booklabel = "12.6.34",
  factor = z^a * exp( -z ),
  begin = [[1, 1 - a + z]],
  general = [ [( 1 - m ) * (m - 1 - a), ( ( 2 * m ) - 1) - a + z] ],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "J-fraction"
):

create( 'contfrac',
  label = "GA.uincgamma.jfrac.02",
  booklabel = "12.6.35",
  factor = z^(a-1) * exp( -z ),
  front = 1,
  begin = [[a - 1, 2 + z - a]],
  general = [[(1 - m)*(m - a), 2*m + z - a]],
  parameters = { a },
  constraints = { abs(functions:-argument(z)) < Pi },
  function = { functions:-uincgamma, GAMMA },
  lhs = functions:-uincgamma(a, z),
  category = "J-fraction"
):

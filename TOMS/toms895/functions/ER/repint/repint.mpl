create( 'series',
  label = "ER.erfcI.asymp.01",
  booklabel = "13.3.2",
  factor = 2 / sqrt(Pi) * exp(-z^2) / (2*z)^( alpha + 1 ),
  general = [( -1 )^k * ( 2*k + alpha )! / ( alpha! * k! * ( 2*z )^( 2*k ) )],
  parameters = { alpha },
  constraints = { abs(functions:-argument(z)) < 3*Pi/4 },
  function = { functions:-erfcI, erfc },
  lhs = functions:-erfcI(alpha, z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "ER.erfcI.sfrac.01",
  booklabel = "13.3.5",
  begin = [[1 / 2, z]],
  general = [[(alpha + m - 1) / 2, z]],
  parameters = { alpha },
  constraints = { Re(z) > 0, alpha::nonnegint },
  function = { functions:-erfcI, erfc },
  lhs = functions:-erfcI(alpha, z) / functions:-erfcI(alpha-1, z),
  category = "S-fraction"
):

create( 'series',
  label = "ER.erfc.power.01",
  booklabel = "13.2.9",
  factor = exp( -z^2 ),
  general = [( -z )^k / GAMMA(k/2 + 1)],
  function = erfc,
  lhs = erfc(z),
  category = "power series"
):

create( 'series',
  label = "ER.cerf.power.01",
  booklabel = "13.2.10",
  general = [( I*z )^k / GAMMA(k/2 + 1)],
  function = { functions:-cerf, erfc },
  lhs = functions:-cerf(z),
  category = "power series"
):

create( 'series',
  label = "ER.erfc.asymp.01",
  booklabel = "13.2.11",
  factor = 1 / ( sqrt(Pi) * z * exp( z^2 ) ),
  general = [( -1 )^k * pochhammer(1 / 2, k) * z^(-2*k)],
  constraints = { abs(functions:-argument(z)) < 3*Pi/4 },
  function = erfc,
  lhs = erfc(z),
  category = "asymptotic series"
):

create( 'series',
  label = "ER.cerf.asymp.01",
  booklabel = "13.2.12",
  factor = I / ( Pi * z ),
  general = [ pochhammer(1 / 2, k) * z^(-2*k)],
  constraints = { abs(functions:-argument(-I*z)) < 3*Pi/4 },
  function =  { functions:-cerf, erfc },
  lhs = functions:-cerf(z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "ER.erfc.sfrac.01",
  booklabel = "13.2.20a",
  factor = z / sqrt(Pi) * exp(-z^2),
  begin = [[1, z^2]],
  even = [(m - 1 ) / 2, 1 ],
  odd = [(m - 1 ) / 2, z^2 ],
  constraints = { Re(z) > 0 },
  function = erfc,
  lhs = erfc(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "ER.cerf.sfrac.01",
  booklabel = "13.2.20b",
  factor = -I*z / sqrt(Pi),
  begin = [[1, -z^2]],
  even = [(m - 1 ) / 2, 1 ],
  odd = [(m - 1 ) / 2, -z^2 ],
  constraints = { Im(z) > 0 },
  function = { functions:-cerf, erfc },
  lhs = functions:-cerf(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "ER.erfc.jfrac.01",
  booklabel = "13.2.23a",
  factor = exp(-z^2) / sqrt(Pi),
  begin = [[2*z, 2*z^2+1]],
  general = [[-(2*m - 3) * (2*m - 2), 4*m - 3 + 2*z^2]],
  constraints = { Re(z) > 0 },
  function = erfc,
  lhs = erfc(z),
  category = "J-fraction"
):

create( 'contfrac',
  label = "ER.cerf.jfrac.01",
  booklabel = "13.2.23b",
  begin = [[-I*z / sqrt(Pi), 1/2 - z^2]],
  general = [[-(- 3/2 + m) * (m - 1), 2 * m - 3/2 - z ]],
  constraints = { Im(z) > 0 },
  function = { functions:-cerf, erfc },
  lhs = functions:-cerf(z),
  category = "J-fraction"
):


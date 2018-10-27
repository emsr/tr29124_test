create( 'series',
  label = "ER.erf.power.01",
  booklabel = "13.1.7",
  factor = 2 / sqrt(Pi),
  general = [( ( -1 )^k * z^( 2*k + 1 ) ) / ( (2*k + 1 ) * k! )],
  function = erf,
  lhs = erf(z),
  category = "power series"
):

create( 'series',
  label = "ER.erf.power.02",
  booklabel = "13.1.8",
  factor = 2 / sqrt(Pi) * exp(-z^2),
  general = [( z^( 2*k + 1 ) ) / ( pochhammer(3/2,k) )],
  function = erf,
  lhs = erf(z),
  category = "power series"
):

create( 'series',
  label = "ER.dawson.power.01",
  booklabel = "13.1.9",
  factor = exp(-z^2),
  general = [z^( 2*k + 1 ) / ( (2*k + 1 ) * k! )],
  function = dawson,
  lhs = dawson(z),
  category = "power series"
):

create( 'series',
  label = "ER.dawson.power.02",
  booklabel = "13.1.10",
  factor = 1,
  general = [( (-1)^k * z^( 2*k + 1 ) ) / ( pochhammer(3/2,k) )],
  function = dawson,
  lhs = dawson(z),
  category = "power series"
):

create( 'contfrac',
  label = "ER.erf.cfrac.01",
  booklabel = "13.1.11a",
  factor = 1 / ( sqrt(Pi) * z * exp(z^2) ),
  begin = [[2 * z^2, 1]],
  odd = [( 2 *( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  even = [( -2 * ( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  function = erf,
  lhs = erf(z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "ER.dawson.cfrac.01",
  booklabel = "13.1.11b",
  factor = -1/(2*z),
  begin = [[-2 * z^2, 1]],
  odd = [ -( 2 *( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  even = [( 2 * ( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  function = dawson,
  lhs = dawson(z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "ER.erf.tfrac.01",
  booklabel = "13.1.13a",
  factor = 1 / ( sqrt(Pi) * z * exp(z^2) ),
  begin = [[2 * z^2, 1 - 2 * z^2]],
  general = [[( 4 * ( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1 -2 / ( 2*m -1 ) * z^2]],
  function = erf,
  lhs = erf(z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "ER.dawson.tfrac.01",
  booklabel = "13.1.13b",
  factor = -1 / (2 * z),
  begin = [[-2 * z^2, 1 + 2 * z^2]],
  general = [[-( 4 * ( m - 1 ) ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1 + 2 / ( 2*m -1 ) * z^2]],
  function = dawson,
  lhs = dawson(z),
  category = "T-fraction"
):


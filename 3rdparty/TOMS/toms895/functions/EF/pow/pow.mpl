create( 'contfrac',
  label = "EF.pow.cfrac.01",
  booklabel = "11.7.1",
  front = 1,
  begin = [[alpha*z, 1]],
  even = [( ( m / 2 ) - alpha ) * z / (2 * (m - 1) ), 1],
  odd = [( ( ( m - 1 ) / 2 ) + alpha ) * z / ( 2 * m ), 1],
  parameters = { alpha },
  constraints = { abs(functions:-argument(z+1)) < Pi },
#  function = pow,
  function = `^`,
  lhs = (1+z)^alpha,
  category = "C-fraction"
):

create( 'contfrac',
  label = "EF.pow.cfrac.02",
  booklabel = "11.7.2",
  begin = [[1, 1], [-alpha * z, 1]],
  general = [ [ ((m-1)/2 + alpha) * z / (2*(m-2)), 1], [(m/2 - 1 - alpha) * z/ (2*(m-1)) , 1] ],
  parameters = { alpha },
  constraints = { abs(functions:-argument(z+1)) < Pi },
#  function = pow,
  function = `^`,
  lhs = (1+z)^alpha,
  category = "C-fraction"
):

create( 'contfrac',
  label = "EF.pow.cfrac.03",
  booklabel = "11.7.3",
  begin = [[1, 1], [-alpha*z, 1+z], [(alpha-1)*z/2, 1]],
  general = [[(-alpha - ((m-2)/2))*z / (2 * (m - 1) * (1+z) ), 1 ], [ (alpha - ((m-1)/2)) * z / (2 * (m - 2) * (1 + z)), 1]],
  parameters = { alpha },
  constraints = { abs(functions:-argument(z+1)) < Pi },
#  function = pow,
  function = `^`,
  lhs = (1+z)^alpha,
  category = "C-fraction"
):

create( 'contfrac',
  label = "EF.pow.cfrac.04",
  booklabel = "11.7.4",
  front = 1,
  begin = [[2 * alpha / z, 1 - alpha / z]],
  general = [ [(alpha^2 - (m-1)^2 ) / ( ( 2 * (m - 1) - 1 ) * (2 * ( m - 1) + 1 ) * z^2), 1] ],
  parameters = { alpha },
  constraints = { ( Re(z) < -1 ) or ( Re(z) > 1 ) },
#  function = pow,
  function = `^`,
  lhs = ((z+1)/(z-1))^alpha,
  category = "C-fraction"
):


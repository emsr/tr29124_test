create( 'contfrac',
  label = "CN.e.rfrac.03",
  booklabel = "10.4.1",
  front = 1,
  begin = [[1, 1]],
  general = [[1, 1],[1, 1],[1, 4 * ( m - 1) / 3 + 1]],
  function = exp,
  lhs = sqrt(exp(1)),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.e.rfrac.04",
  booklabel = "10.4.2",
  front = 1,
  general = [[1, ( 2 * ( ( m + 2 ) / 3 ) - 1 ) * alpha - 1], [1, 1], [1, 1]],
  parameters = { alpha },
  function = exp,
  lhs = exp(1/alpha),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.e.rfrac.05",
  booklabel = "10.4.3",
  front = (alpha + 1) / alpha,
  factor = 1/alpha,
  general = [ [1, 2 * alpha - 1], [1, (2/3) * (m+1) ], [1, 1] ],
  parameters = { alpha },
  function = exp,
  lhs = exp(1/alpha),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.e.rfrac.06",
  booklabel = "10.4.4",
  factor = alpha,
  begin = [ [1, alpha-1], [1, 2*alpha] ],
  general = [ [1, 1], [1, (2/3) * (m-1) ], [1, 2*alpha - 1] ],
  parameters = { alpha },
  function = exp,
  lhs = exp(1/alpha),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.e.rfrac.07",
  booklabel = "10.4.5",
  front = 7,
  general = [ [ 1, 3 * (m+4)/5 - 1 ], [ 1, 1 ], [ 1, 1 ], [ 1, 3 * (m+1)/5 ], [ 1, 12*m/5 + 6 ] ],
  function = exp,
  lhs = exp(1)*exp(1),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.pi.sfrac.03",
  booklabel = "10.4.6",
  begin = [[1, 1]],
  general = [ [(m-1)^4, 2 * m - 1] ],
  function = Pi,
  lhs = Pi^2 / 12,
  category = "S-fraction"
):

create( 'contfrac',
  label = "CN.e.sfrac.01",
  booklabel = "10.4.7",
  begin = [[alpha, beta]],
  general = [ [alpha^2, ( 2 * m - 1 ) * beta] ],
  parameters = { alpha, beta },
  function = exp,
  lhs = ( exp(2*alpha/beta) - 1 ) / ( exp(2*alpha/beta) + 1 ),
  category = "S-fraction"
):

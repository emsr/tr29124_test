create( 'series',
  label = "CN.e.power.01",
  booklabel = "10.3.1b",
  general = [ 1 / ( k)! ],
  function = exp,
  lhs = exp(1),
  category = "power series"
):

create( 'contfrac',
  label = "CN.e.rfrac.01",
  booklabel = "10.3.5",
  front = 2,
  general = [[1, 1], [1 , 2 / 3 * ( m + 1 )], [1, 1]],
  function = exp,
  lhs = exp(1),
  category = "regular continued fraction"
):

create( 'contfrac',
  label = "CN.e.rfrac.02",
  booklabel = "10.3.6",
  general = [ [1, ( 4 * m - 2 )] ],
  function = exp,
  lhs = ( exp(1) - 1 ) / ( exp(1) + 1 ),
  category = "regular continued fraction"
):

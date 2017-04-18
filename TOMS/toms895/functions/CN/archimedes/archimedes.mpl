create( 'series',
  label = "CN.pi.power.01",
  booklabel = "10.2.1",
  factor = 4,
  general = [ ( -1 )^( k ) / ( 2 * k + 1 ) ],
  function = Pi,
  lhs = Pi,
  category = "power series"
):

create( 'contfrac',
  label = "CN.pi.sfrac.01",
  booklabel = "10.2.5",
  begin = [[4, 1]],
  general = [ [( m - 1 )^2, 2 * m - 1] ],
  function = Pi,
  lhs = Pi,
  category = "S-fraction"
):

create( 'contfrac',
  label = "CN.pi.sfrac.02",
  booklabel = "10.2.6",
  front = 3,
  general = [ [( 2 * m - 1 )^2, 6] ],
  function = Pi,
  lhs = Pi,
  category = "S-fraction"
):

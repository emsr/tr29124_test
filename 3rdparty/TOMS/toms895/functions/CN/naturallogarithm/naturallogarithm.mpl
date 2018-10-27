create( 'series',
  label = "CN.ln2.power.01",
  booklabel = "10.5.2",
  general = [ ( ( -1 )^(k+1) * 1^( k ) ) / ( k + 1 ) ],
  function = ln,
  lhs = ln(2),
  category = "power series"
):

create( 'contfrac',
  label = "CN.ln2.sfrac.01",
  booklabel = "10.5.3",
  begin = [[1,1]],
  general = [ [( m - 1 )^2, 1] ],
  function = ln,
  lhs = ln(2),
  category = "S-fraction"
):

create( 'contfrac',
  label = "CN.ln2.sfrac.02",
  booklabel = "10.5.4",
  begin = [[1,1]],
  even = [m / ( 4 * m - 4 ), 1],
  odd = [( m - 1 ) / ( 4 * m ), 1],
  function = ln,
  lhs = ln(2),
  category = "S-fraction"
):

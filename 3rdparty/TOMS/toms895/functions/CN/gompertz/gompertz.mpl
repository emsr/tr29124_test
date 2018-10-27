create( 'contfrac',
  label = "CN.gompertz.frac.01",
  booklabel = "10.13.1",
  begin = [[1, 2]],
  general = [ [-( m - 1 )^2, 2 * m ] ],
  function = functions:-Gompertz,
  lhs = functions:-Gompertz(),
  category = "continued fraction"
):

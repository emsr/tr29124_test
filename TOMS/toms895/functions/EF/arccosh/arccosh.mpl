create( 'series',
  label = "EF.arccosh.power.01",
  booklabel = "11.6.2",
  front = ln(2/z),
  general = [ -(doublefactorial(2*k-1)/(doublefactorial(2*k)*(2*k)))*z^(2*k) ],
  constraints = { abs(z) < 1 },
  function = arccosh,
  lhs = arccosh(1/z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.arccosh.sfrac.01",
  booklabel = "11.6.6",
  begin = [[z*sqrt(z^2-1), 1]],
  even = [(m*(m-1)/((2*m-3)*(2*m-1)))*(z^2-1), 1],
  odd = [((m-1)*(m-2)/((2*m-3)*(2*m-1)))*(z^2-1), 1],
  constraints = { Re(z) > 0 },
  function = arccosh,
  lhs = arccosh(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arccosh.sfrac.02",
  booklabel = "11.6.7",
  begin = [[sqrt(z^2-1)/z, 1]],
  general = [ [-((m-1)^2/((2*m-3)*(2*m-1)))*(z^2-1)/z^2, 1] ],
  constraints = { abs(functions:-argument(1/z^2)) < Pi },
  function = arccosh,
  lhs = arccosh(z),
  category = "S-fraction"
):


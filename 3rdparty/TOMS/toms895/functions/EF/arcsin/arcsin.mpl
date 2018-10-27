create( 'series',
  label = "EF.arcsin.power.01",
  booklabel = "11.4.1",
  front = z,
  general = [ ((doublefactorial(2*k-1))/(doublefactorial(2*k)*(2*k+1)))*z^(2*k+1) ],
  constraints = { abs(z) < 1 },
  function = arcsin,
  lhs = arcsin(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.arcsin.sfrac.01",
  booklabel = "11.4.4",
  begin = [[z/sqrt(1-z^2), 1]],
  general = [ [((m-1)^2/((2*m-3)*(2*m-1)))*z^2/(1-z^2), 1] ],
  constraints = { abs(functions:-argument(1-z^2)) < Pi },
  function = arcsin,
  lhs = arcsin(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arcsin.sfrac.02",
  booklabel = "11.4.5",
  begin = [[z*sqrt(1-z^2), 1]],
  even = [(-m*(m-1)/((2*m-1)*(2*m-3)))*z^2, 1],
  odd = [(-(m-1)*(m-2)/((2*m-1)*(2*m-3)))*z^2, 1],
  constraints = { abs(functions:-argument(1-z^2)) < Pi },
  function = arcsin,
  lhs = arcsin(z),
  category = "S-fraction"
):

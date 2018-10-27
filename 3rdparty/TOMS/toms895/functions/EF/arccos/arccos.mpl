create( 'series',
  label = "EF.arccos.power.01",
  booklabel = "11.4.2",
  front = Pi/2 - z,
  general = [ -((doublefactorial(2*k-1))/(doublefactorial(2*k)*(2*k+1)))*z^(2*k+1) ],
  constraints = { abs(z) < 1 },
  function = arccos,
  lhs = arccos(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.arccos.sfrac.01",
  booklabel = "11.4.6",
  begin = [[sqrt(1-z^2) / z, 1]],
  general = [ [((m-1)^2/((2*m-3)*(2*m-1)))*(1-z^2)/z^2, 1] ],
  constraints = { Re(z) > 0, (1-z^2)/z^2 > -1 },
  function = arccos,
  lhs = arccos(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arccos.sfrac.02",
  booklabel = "11.4.7",
  begin = [[z*sqrt(1-z^2), 1]],
  even = [(-m*(m-1)/((2*m-1)*(2*m-3)))*(1-z^2), 1],
  odd = [(-(m-1)*(m-2)/((2*m-1)*(2*m-3)))*(1-z^2), 1],
  constraints = { Re(z) > 0 },
  function = arccos,
  lhs = arccos(z),
  category = "S-fraction"
):

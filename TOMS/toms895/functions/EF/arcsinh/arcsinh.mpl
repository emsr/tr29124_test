create( 'series',
  label = "EF.arcsinh.power.01",
  booklabel = "11.6.1",
  front = z,
  general = [ (((-1)^k*doublefactorial(2*k-1))/(doublefactorial(2*k)*(2*k+1)))*z^(2*k+1) ],
  constraints = { abs(z) < 1 },
  function = arcsinh,
  lhs = arcsinh(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.arcsinh.sfrac.01",
  booklabel = "11.6.4",
  begin = [[z*sqrt(1+z^2), 1]],
  even = [(m*(m-1)/((2*m-3)*(2*m-1)))*z^2, 1],
  odd = [(((m-1)*(m-2))/((2*m-3)*(2*m-1)))*z^2, 1],
  constraints = { not(I*z < -1), not(I*z > 1) },
  function = arcsinh,
  lhs = arcsinh(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arcsinh.sfrac.02",
  booklabel = "11.6.5",
  begin = [[z/sqrt(1+z^2), 1]],
  general = [ [-(((m-1)^2)/((2*m-3)*(2*m-1)))*z^2/(1+z^2), 1] ],
  constraints = { not(I*z < -1), not(I*z > 1) },
  function = arcsinh,
  lhs = arcsinh(z),
  category = "S-fraction"
):

create( 'series',
  label = "EF.arctan.power.01",
  booklabel = "11.4.3",
  general = [ ((-1)^(k)/(2*k+1))*z^(2*k+1) ],
  constraints = { abs(z) < 1 },
  function = arctan,
  lhs = arctan(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.arctan.sfrac.01",
  booklabel = "11.4.8",
  begin = [[z, 1]],
  general = [ [(m-1)^2*z^2, 2*m-1] ],
  constraints = { not(I*z < -1), not(I*z > 1) },
  function = arctan,
  lhs = arctan(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arctan.sfrac.02",
  booklabel = "11.4.9",
  begin = [[z/(1+z^2), 1]],
  even = [-((m*(m-1))/((2*m-3)*(2*m-1)))*(z^2/(1+z^2)),1],
  odd = [-(((m-1)*(m-2))/((2*m-3)*(2*m-1)))*(z^2/(1+z^2)),1],
  constraints = { not(I*z < -1), not(I*z > 1) },
  function = arctan,
  lhs = arctan(z),
  category = "S-fraction"
):

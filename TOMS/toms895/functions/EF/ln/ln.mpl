create( 'series',
  label = "EF.ln.power.01",
  booklabel = "11.2.1",
  general = [ ((-1)^(k+2)/(k+1))*z^(k+1) ],
  constraints = { abs(z) < 1 },
  function = ln,
  lhs = ln(1+z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.ln.sfrac.01",
  booklabel = "11.2.2",
  begin = [[z, 1]],
  even = [(m/(4*(m-1)))*z, 1],
  odd = [((m-1)/(4*m))*z, 1],
  constraints = { abs(functions:-argument(1+z)) < Pi },
  function = ln,
  lhs = ln(1+z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.ln.e2frac.01",
  booklabel = "11.2.3",
  begin = [[2*z, 2+z]],
  general = [ [-(m-1)^2*z^2,(2*m-1)*(2+z)] ],
  constraints = { abs(functions:-argument(1-z^2/(2+z)^2)) < Pi },
  function = ln,
  lhs = ln(1+z),
  category = "Euler even contraction"
):

create( 'contfrac',
  label = "EF.ln.sfrac.02",
  booklabel = "11.2.4",
  begin = [[2*z, 1]],
  general = [ [(-((m-1)^2/((2*m-3)*(2*m-1))) * z^2),1] ],
  constraints = { abs(functions:-argument(1-z^2)) < Pi },
  function = ln,
  lhs = ln((1+z)/(1-z)),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.ln.thiele.01",
  booklabel = "6.8.8",
  front = z,
  begin = [[-z^2/2, 1]],
  even = [((m/2) + 1)^2 * z / (m*(m+1)), 1],
  odd = [((m-1)/2)^2 * z / (m*(m+1)), 1],
  function = ln,
  lhs = ln(1+z),
  category = "Thiele interpolating fraction"
):

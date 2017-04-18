create( 'series',
  label = "EF.coth.power.01",
  booklabel = "11.5.4",
  general = [ (4^k*bernoulli(2*k)/(2*k)!)*z^(2*k-1) ],
  constraints = { abs(z) < Pi },
  function = coth,
  lhs = coth(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.coth.thiele.01",
  booklabel = "11.5.6",
  front = 1/z,
  begin = [[4*Pi^(-2)*z, 1]],
  general = [ [(m-1)^2*((m-1)^2+4*Pi^(-2)*z^2), 2*m-1] ],
  function = coth,
  lhs = coth(z),
  category = "Thiele interpolating fraction"
):

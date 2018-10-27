create( 'series',
  label = "EF.tanh.power.01",
  booklabel = "11.5.3",
  general = [ (4^(k+1)*(4^(k+1)-1)*(bernoulli(2*(k+1))/(2*(k+1))!))*z^(2*(k+1)-1) ],
  constraints = { abs(z) < Pi/2 },
  function = tanh,
  lhs = tanh(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.tanh.sfrac.01",
  booklabel = "11.5.5",
  begin = [[z, 1]],
  general = [ [(1/((2*m-3)*(2*m-1)))*z^2, 1] ],
  function = tanh,
  lhs = tanh(z),
  category = "S-fraction"
):

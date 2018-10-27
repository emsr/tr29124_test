create( 'series',
  label = "EF.exp.power.01",
  booklabel = "11.1.1",
  general = [ z^k/k! ],
  function = exp,
  lhs = exp(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.exp.sfrac.01",
  booklabel = "11.1.2",
  front = 1,
  begin = [[2*z, 2-z], [z^2/6, 1]],
  general = [ [(1/(4*(2*m-3)*(2*m-1)))*z^2, 1] ],
  function = exp,
  lhs = exp(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.exp.cfrac.01",
  booklabel = "11.1.3",
  front = 1,
  begin = [[z, 1]],
  even = [(-1/(2*(m-1)))*z, 1],
  odd = [(1/(2*m))*z, 1],
  function = exp,
  lhs = exp(z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "EF.exp.tfrac.01",
  booklabel = "11.1.4",
  front = 1,
  begin = [[z,1-z]],
  general = [ [(m-1)*z, m-z] ],
  function = exp,
  lhs = exp(z),
  category = "T-fraction"
):

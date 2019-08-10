create( 'series',
  label = "CN.zeta.power.01",
  booklabel = "10.11.1",
  general = [ 1 / ( k + 1 )^z ],
  function = Zeta,
  lhs = Zeta(z),
  category = "power series"
):

create( 'series',
  label = "CN.zeta3.power.01",
  booklabel = "10.11.3",
  general = [ ( -1 )^( k) * ( ( (k!)^10 * ( 205 * k^2 + 250 * k + 77) ) / ( 64 * ( ( 2*k + 1 )! )^5 ) ) ],
  function = Zeta,
  lhs = Zeta(3),
  category = "power series"
):

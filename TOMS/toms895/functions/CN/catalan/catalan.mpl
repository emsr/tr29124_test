create( 'series',
  label = "CN.catalan.power.01",
  booklabel = "10.12.1",
  general = [ (-1)^( k ) / ( 2 * k + 1 )^2  ],
  function = Catalan,
  lhs = Catalan,
  category = "power series"
):

create( 'contfrac',
  label = "CN.catalan.sfrac.01",
  booklabel = "10.12.3",
  factor = 1 / 2,
  front = 1 / 2,
  begin = [[1 , 1/2]],
  even = [( m / 2 )^2 , 1/2],
  odd = [ ( ( m - 1 ) / 2  ) * ( ( m + 1 ) / 2 ) , 1/2],
  function = Catalan,
  lhs = Catalan,
  category = "S-fraction"
):

create( 'contfrac',
  label = "CN.catalan.frac.01",
  booklabel = "10.12.5",
  begin = [[13/2, 7]],
  general = [[(2*m-3)^4 * (2*m-2)^4 * (20*(m-2)^2-8*(m-2)+1)* (20*(m)^2-8*(m)+1), 3520*(m-1)^6+5632*(m-1)^5+2064*(m-1)^4-384*(m-1)^3-156*(m-1)^2+16*(m-1)+7]],
  function = Catalan,
  lhs = Catalan,
  category = "continued fraction"
):


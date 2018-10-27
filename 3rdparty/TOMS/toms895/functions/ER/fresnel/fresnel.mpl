create( 'series',
  label = "ER.fresnelC.power.01",
  booklabel = "13.4.6a",
  general = [( ( -1 )^k * ( Pi/2 )^( 2*k ) ) / ( ( 2*k )! * ( 4*k + 1 ) ) * z^( 4*k + 1 )],
  function = FresnelC,
  lhs = FresnelC(z),
  category = "power series"
):

create( 'series',
  label = "ER.fresnelS.power.01",
  booklabel = "13.4.6b",
  general = [( ( -1 )^k * ( Pi/2 )^( 2*k + 1 ) ) / ( ( 2*k + 1 )! * ( 4*k + 3 ) ) * z^( 4*k + 3 )],
  function = FresnelS,
  lhs = FresnelS(z),
  category = "power series"
):

create( 'contfrac',
  label = "ER.fresnel.cfrac.01",
  booklabel = "13.4.9",
  factor = exp( ( I * Pi * z^2 ) / 2 ) / z,
  begin = [[z^2, 1]],
  even = [I * Pi * ( m - 1 ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  odd = [-I * Pi * (  m - 1 ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1],
  function = { FresnelC, FresnelS },
  lhs = FresnelC(z) + I*FresnelS(z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "ER.fresnel.tfrac.01",
  booklabel = "13.4.10",
  factor = exp( ( I * Pi * z^2 ) / 2 ) / z,
  begin = [[z^2, 1 + I * Pi * z^2]],
  general = [[-2 * I * Pi * ( m - 1 ) / ( ( 2*m - 3 ) * ( 2*m - 1 ) ) * z^2, 1 + I * Pi / ( 2*m - 1 ) * z^2 ]],
  function = { FresnelC, FresnelS },
  lhs = FresnelC(z) + I*FresnelS(z),
  category = "T-fraction"
):

create( 'series',
  label = "BS.besselI.power.01",
  booklabel = "17.2.20",
  factor = (z/2)^nu,
  general = [  (z/2)^(2*k) / (k! * GAMMA(nu+k+1)) ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselI,
  lhs = BesselI(nu, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselI.power.02",
  booklabel = "17.2.21",
  factor = exp(-z) * (z/2)^nu / GAMMA(nu+1),
  general = [ pochhammer(nu+1/2, k) * 2^k * z^k / (pochhammer(2*nu+1, k) * k!)  ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselI,
  lhs = BesselI(nu, z),
  category = "power series"
):

create( 'series',
  label = "BS.besseli.power.01",
  booklabel = "17.2.22",
  factor = sqrt(Pi) / ((2*n+1) * GAMMA(n+1/2)) * (z/2)^n,
  general = [ ( z^2/4 )^k / (k! * pochhammer(n+3/2, k)) ],
  parameters = {n},
  constraints = { n::integer },
  function = { functions:-Besseli, BesselI },
  lhs = functions:-Besseli(n, z),
  category = "power series"
):

create( 'series',
  label = "BS.besseli.power.02",
  booklabel = "17.2.23",
  factor = sqrt(Pi) * exp(-I*z) / ((2*n+1) * GAMMA(n+1/2)) * (z/2)^n,
  general = [ pochhammer(n+1, k) * ( 2*z )^k / (k! * pochhammer(2*n+2, k)) ],
  parameters = {n},
  constraints = { n::integer },
  function = { functions:-Besseli, BesselI },
  lhs = functions:-Besseli(n, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselI.asymp.01",
  booklabel = "17.2.24",
  general = [ ((-1)^k * exp(z) + exp(-z + (2*nu + 1) * I * Pi/2)) * functions:-hankelsymb(nu, k) / (sqrt(2*Pi*z) * (2*z)^k) ],
  parameters = {nu},
  constraints = { functions:-argument(z)::RealRange(Open(-Pi/2),Open(3*Pi/2)) },
  function = BesselI,
  lhs = BesselI(nu, z),
  category = "asymptotic series"
):

create( 'series',
  label = "BS.besselI.asymp.02",
  booklabel = "17.2.25",
  general = [ ((-1)^k * exp(z) + exp(-z - (2*nu + 1) * I * Pi/2)) * functions:-hankelsymb(nu, k) / (sqrt(2*Pi*z) * (2*z)^k) ],
  parameters = {nu},
  constraints = { functions:-argument(z)::RealRange(Open(-3*Pi/2),Open(Pi/2)) },
  function = BesselI,
  lhs = BesselI(nu, z),
  category = "asymptotic series"
):

create( 'series',
  label = "BS.besselK.asymp.01",
  booklabel = "17.2.27",
  factor = sqrt(Pi/(2*z)) * exp(-z),
  general = [ functions:-hankelsymb(nu, k) / (-2*z)^k ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < 3*Pi/2 },
  function = BesselK,
  lhs = BesselK(nu, z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "BS.besselI.sfrac.01",
  booklabel = "17.2.32",
  begin = [[z / (2*(nu+1)), 1]],
  general = [[ 1/(4*(nu+m-1)*(nu+m)) * z^2, 1]],
  parameters = {nu},
  constraints = { nu::nonnegative },
  function = BesselI,
  lhs = BesselI(nu+1, z) / BesselI(nu, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "BS.besseli.sfrac.01",
  booklabel = "17.2.33",
  begin = [[z / (2*n+3), 1]],
  general = [[ 1/(4*((n+1/2)+m-1)*((n+1/2)+m)) * z^2, 1]],
  parameters = {n},
  constraints = { n::nonnegint },
  function = { functions:-Besseli, BesselI },
  lhs = functions:-Besseli(n+1, z) / functions:-Besseli(n, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "BS.besselK.cfrac.01",
  booklabel = "17.2.34",
  front = nu/z,
  factor = -1,
  begin = [[1, 1], [(-2*nu-1)/(2*z), 1]],
  general = [[ ( m/2 + nu ) / (2*z), 1 ], [ ( (m-3)/2 - nu ) / (2*z), 1 ]],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselK,
  lhs = nu/z - BesselK(nu+1, z) / BesselK(nu, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "BS.besselI.tfrac.01",
  booklabel = "17.2.38",
  begin = [[z, 2*nu + 2 + z]],
  general = [[ -(2*nu + 2*m - 1) * z, 2*nu+m+1 + 2*z]],
  parameters = {nu},
  constraints = { nu::Not(negint) },
  function = BesselI,
  lhs = BesselI(nu+1, z) / BesselI(nu, z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "BS.besseli.tfrac.01",
  booklabel = "17.2.39",
  begin = [[z, 2*n + 3 + z]],
  general = [[ -2*(n + m) * z, 2*n+m+2 + 2*z]],
  parameters = {n},
  constraints = { n::nonnegint },
  function = { functions:-Besseli, BesselI },
  lhs = functions:-Besseli(n+1, z) / functions:-Besseli(n, z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "BS.besselK.jfrac.01",
  booklabel = "17.2.40",
  front = nu/z-(2*nu+1+2*z) / (2*z),
  factor = -1/z,
  general = [[ nu^2 - (2*m-1)^2 / 4, 2 * (z+ m) ]],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselK,
  lhs = nu/z - BesselK(nu+1, z) / BesselK(nu, z),
  category = "J-fraction"
):

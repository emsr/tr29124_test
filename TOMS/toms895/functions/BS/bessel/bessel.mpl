create( 'series',
  label = "BS.besselJ.power.01",
  booklabel = "17.1.2a",
  factor = (z/2)^nu,
  general = [ (-1)^k / ( k! * GAMMA(nu+k+1) ) * (z/2)^(2*k) ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselJ,
  lhs = BesselJ(nu, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselj.power.01",
  booklabel = "17.1.12a",
  factor = sqrt(Pi/(2*z)) * (z/2)^(n+1/2),
  general = [ (-1)^k / ( k! * GAMMA((n+1/2)+k+1) ) * (z/2)^(2*k) ],
  parameters = {n},
  constraints = { n::integer },
  function = { functions:-Besselj, BesselJ },
  lhs = functions:-Besselj(n, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselJ.power.02",
  booklabel = "17.1.22",
  factor = exp(-I*z) * (z/2)^nu / GAMMA(nu+1),
  general = [pochhammer(nu+1/2, k) * (2*I)^k *z^k / (pochhammer(2*nu+1, k) * k!)],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselJ,
  lhs = BesselJ(nu, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselj.power.02",
  booklabel = "17.1.25",
  factor = sqrt(Pi) / ((2*n+1) * GAMMA(n+1/2)) * (z/2)^n,
  general = [1/pochhammer(n+3/2, k) * (-z^2/4)^k / k! ],
  parameters = {n},
  constraints = { n::integer },
  function = { functions:-Besselj, BesselJ },
  lhs = functions:-Besselj(n, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselj.power.03",
  booklabel = "17.1.26",
  factor = sqrt(Pi) * exp(-I*z) / ((2*n+1) * GAMMA(n+1/2)) * (z/2)^n,
  general = [pochhammer(n+1, k)/pochhammer(2*n+2, k) * (2*I*z)^k / k! ],
  parameters = {n},
  constraints = { n::integer },
  function = { functions:-Besselj, BesselJ },
  lhs = functions:-Besselj(n, z),
  category = "power series"
):

create( 'series',
  label = "BS.besselJ.asymp.01",
  booklabel = "17.1.28",
  factor = sqrt(2/(Pi*z)),
  general = [ (-1)^k * functions:-hankelsymb(nu, 2*k) / (2*z)^(2*k) * cos(z - (nu/2 + 1/4) * Pi) - (-1)^k * functions:-hankelsymb(nu, 2*k+1) / (2*z)^(2*k+1) * sin(z - (nu/2 + 1/4) * Pi) ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselJ,
  lhs = BesselJ(nu, z),
  category = "asymptotic series"
):

create( 'series',
  label = "BS.besselY.asymp.01",
  booklabel = "17.1.29",
  factor = sqrt(2/(Pi*z)),
  general = [ (-1)^k * functions:-hankelsymb(nu, 2*k) / (2*z)^(2*k) * sin(z - (nu/2 + 1/4) * Pi) + (-1)^k * functions:-hankelsymb(nu, 2*k+1) / (2*z)^(2*k+1) * cos(z - (nu/2 + 1/4) * Pi) ],
  parameters = {nu},
  constraints = { abs(functions:-argument(z)) < Pi },
  function = BesselY,
  lhs = BesselY(nu, z),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "BS.besselJ.sfrac.01",
  booklabel = "17.1.38",
  begin = [[z / (2*nu + 2), 1]],
  general = [[(I*z)^2 / (4*(nu+m-1)*(nu+m)), 1]],
  parameters = {nu},
  constraints = { nu::nonnegative },
  function = BesselJ,
  lhs = BesselJ(nu+1, z) / BesselJ(nu, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "BS.besselj.sfrac.01",
  booklabel = "17.1.39",
  begin = [[z/(2*n+3), 1]],
  general = [[(I*z)^2 / (4 * (n+1/2 + m - 1) * (n+1/2 + m) ), 1]],
  parameters = {n},
  constraints = { n::nonnegint },
  function = { functions:-Besselj, BesselJ },
  lhs = functions:-Besselj(n+1, z) / functions:-Besselj(n, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "BS.besselJ.sfrac.02",
  booklabel = "17.1.40",
  factor = -1,
  general = [[-1, 2*(nu+m)/z]],
  parameters = {nu},
  constraints = { nu::nonnegint },
  function = BesselJ,
  lhs = BesselJ(nu+1, z) / BesselJ(nu, z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "BS.hankelH1.cfrac.01",
  booklabel = "17.1.44",
  begin = [[ -1, 1 ]],
  even = [ (m-3-2*nu) / (-2*I*z), 1 ],
  odd = [ (m+2*nu) / (-2*I*z), 1],
  parameters = {nu},
  constraints = { abs(functions:-argument(I*z)) < Pi },
  function = HankelH1,
  lhs = HankelH1(nu+1, z) / HankelH1(nu, z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "BS.besselJ.tfrac.01",
  booklabel = "17.1.48",
  begin = [[z, 2*nu + 2 - I*z]],
  general = [[(2*nu + 2*m - 1)*I * z, 2*nu + m + 1 + (-2*I) *z]],
  parameters = {nu},
  constraints = { nu::Not(negint) },
  function = BesselJ,
  lhs = BesselJ(nu+1, z) / BesselJ(nu, z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "BS.besselj.tfrac.01",
  booklabel = "17.1.49",
  begin = [[z, 2*n + 3 - I*z]],
  general = [[2*(n + m) * I * z, 2*n + m + 2 + (-2*I) *z]],
  parameters = {n},
  constraints = { n::nonnegint },
  function = { functions:-Besselj, BesselJ },
  lhs = functions:-Besselj(n+1, z) / functions:-Besselj(n, z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "BS.hankelH1.jfrac.01",
  booklabel = "17.1.51",
  front = (2*nu + 1 - 2*I*z) / (2*z),
  factor = -1/z,
  general = [ [nu^2 - (2*m-1)^2 / 4, 2*(I*z-m)] ],
  constraints = { abs(functions:-argument(-I*z)) < Pi },
  parameters = {nu},
  function = HankelH1,
  lhs = HankelH1(nu+1, z) / HankelH1(nu, z),
  category = "J-fraction"
):

create( 'series',
  label = "SM.gammacdf.power.01",
  booklabel = "18.4.15",
  factor = (x/theta)^alpha / GAMMA(alpha),
  general = [(-x/theta)^k / ((alpha+k) * k!)],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-gammacdf, # P
  lhs = functions:-gammacdf(x, alpha, theta),
  category = "power series"
):

create( 'series',
  label = "SM.gammacdf.power.02",
  booklabel = "18.4.16",
  factor = (x/theta)^alpha * exp(-x/theta) / GAMMA(alpha+1),
  general = [(x/theta)^k / pochhammer(alpha+1, k)],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-gammacdf, # P
  lhs = functions:-gammacdf(x, alpha, theta),
  category = "power series"
):

create( 'series',
  label = "SM.cgammacdf.asymp.01",
  booklabel = "18.4.17",
  factor = (x/theta)^(alpha-1) * exp(-x/theta) / GAMMA(alpha),
  general = [(-1)^k * pochhammer(1-alpha, k) * (theta/x)^k],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::positive, alpha::RealRange(Open(0), Open(1)), theta::positive },
  function = functions:-cgammacdf, # Q
  lhs = functions:-cgammacdf(x,alpha,theta),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "SM.cgammacdf.sfrac.01",
  booklabel = "18.4.19",
  factor = (x/theta)^alpha * exp(-x/theta) / GAMMA(alpha),
  begin = [[1, x/theta]],
  general = [[m/2-alpha, 1], [(m-1)/2, x/theta]],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::RealRange(Open(0), Open(1)), theta::positive },
  function = functions:-cgammacdf, # Q
  lhs = functions:-cgammacdf(x,alpha,theta),
  category = "S-fraction"
):

create( 'contfrac',
  label = "SM.gammacdf.cfrac.01",
  booklabel = "18.4.20",
  factor = (x/theta)^(alpha-1) * exp(-x/theta) / GAMMA(alpha),
  begin = [[1/alpha * (x/theta), 1]],
  general = [[-(alpha+m/2-1)/ ((alpha+m-2)*(alpha+m-1)) * (x/theta), 1], [((m-1)/2)/ ((alpha+m-2)*(alpha+m-1)) * (x/theta), 1]],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-gammacdf, # P
  lhs = functions:-gammacdf(x,alpha,theta),
  category = "C-fraction"
):

create( 'contfrac',
  label = "SM.cgammacdf.cfrac.01",
  booklabel = "18.4.21",
  front = 1,
  factor = -(x/theta)^(alpha-1) * exp(-x/theta) / GAMMA(alpha),
  begin = [[1/alpha * (x/theta), 1]],
  general = [[-(alpha+m/2-1)/ ((alpha+m-2)*(alpha+m-1)) * (x/theta), 1], [((m-1)/2)/ ((alpha+m-2)*(alpha+m-1)) * (x/theta), 1]],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-cgammacdf, # Q
  lhs = functions:-cgammacdf(x,alpha,theta),
  category = "C-fraction"
):

create( 'contfrac',
  label = "SM.gammacdf.mfrac.01",
  booklabel = "18.4.24",
  begin = [[(x/theta)^alpha * exp(-x/theta) / GAMMA(alpha), alpha - x/theta]],
  general = [[(m-1) * x/theta, alpha+m-1-x/theta ]],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-gammacdf, # P
  lhs = functions:-gammacdf(x,alpha,theta),
  category = "M-fraction"
):

create( 'contfrac',
  label = "SM.cgammacdf.mfrac.01",
  booklabel = "18.4.25",
  front = 1,
  begin = [[-(x/theta)^alpha * exp(-x/theta) / GAMMA(alpha), alpha - x/theta]],
  general = [[(m-1) * x/theta, alpha+m-1-x/theta ]],
  parameters = {alpha, theta},
  variable = x,
  constraints = { x::nonnegative, alpha::positive, theta::positive },
  function = functions:-cgammacdf, # Q
  lhs = functions:-cgammacdf(x,alpha,theta),
  category = "M-fraction"
):

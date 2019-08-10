create( 'series',
  label = "SM.normalcdf.power.01",
  booklabel = "18.2.13a",
  front = 1/2,
  factor = 1/sqrt(2*Pi),
  general = [(-1)^(k-1) * x^(2*k-1) / (2^(k-1) * (2*k-1) * (k-1)!)],
  variable = x,
  constraints = { x::realcons },
  function = functions:-normaldcdf, # F
  lhs = functions:-normaldcdf(x),
  category = "power series"
):

create( 'series',
  label = "SM.cnormalcdf.power.01",
  booklabel = "18.2.13b",
  front = 1/2,
  factor = -1/sqrt(2*Pi),
  general = [(-1)^(k-1) * x^(2*k-1) / (2^(k-1) * (2*k-1) * (k-1)!)],
  variable = x,
  constraints = { x::realcons },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "power series"
):

create( 'series',
  label = "SM.cnormalcdf.asymp.01",
  booklabel = "18.2.14",
  factor = exp(-x^2/2) / (x * sqrt(2*Pi)),
  general = [(-1)^k * doublefactorial(2*k-1)/ x^(2*k)],
  variable = x,
  constraints = { x::realcons },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "SM.cnormalcdf.sfrac.01",
  booklabel = "18.2.16",
  factor = x * exp(-x^2/2) / sqrt(2*Pi),
  begin = [[1, x^2]],
  even = [m-1, 1],
  odd = [m-1, x^2],
  variable = x,
  constraints = { x::positive },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "S-fraction"
):

create( 'contfrac',
  label = "SM.cnormalcdf.sfrac.02",
  booklabel = "18.2.17",
  factor = exp(-x^2/2) / sqrt(2*Pi),
  begin = [[1, x]],
  general = [[m-1, x]],
  variable = x,
  constraints = { x::positive },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "S-fraction"
):

create( 'contfrac',
  label = "SM.cnormalcdf.cfrac.01",
  booklabel = "18.2.18",
  front = 1/2,
  factor = -exp(-x^2/2) / (x*sqrt(2*Pi)),
  begin = [[x^2, 1]],
  general = [[(-1)^(m-1) * (m-1) * x^2 / ((2*m-3)*(2*m-1)), 1]],
  variable = x,
  constraints = { x::realcons },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "C-fraction"
):

create( 'contfrac',
  label = "SM.cnormalcdf.jfrac.01",
  booklabel = "18.2.19",
  factor = exp(-x^2/2) / (sqrt(2*Pi)),
  begin = [[x, 1+x^2]],
  general = [[-(2*m-3)*(2*m-2), 4*m-3+x^2]],
  variable = x,
  constraints = { x::positive },
  function = functions:-cnormaldcdf, # Q
  lhs = functions:-cnormaldcdf(x),
  category = "J-fraction"
):

create( 'contfrac',
  label = "SM.mills.sfrac.01",
  booklabel = "18.2.21",
  begin = [[1, x]],
  general = [[m-1, x]],
  variable = x,
  constraints = { x::positive },
  function = { functions:-Mills, functions:-cnormaldcdf },
  lhs = functions:-Mills(x),
  category = "S-fraction"
):

create( 'contfrac',
  label = "SM.mills.cfrac.01",
  booklabel = "18.2.22",
  front = sqrt(Pi/2) * exp(x^2/2),
  begin = [[-x, 1]],
  general = [[(-1)^(m-1) * (m-1) * x^2, 2*m-1]],
  variable = x,
  constraints = { x::positive },
  function = { functions:-Mills, functions:-cnormaldcdf },
  lhs = functions:-Mills(x),
  category = "C-fraction"
):

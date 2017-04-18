create( 'series',
  label = "SM.cnormalI.asymp.01",
  booklabel = "18.3.3",
  factor = exp(-x^2/2) / (x^(l+1) * sqrt(2*Pi)),
  general = [(-1)^k * (2*k+l)!/ (l! * k! * 2^k * x^(2*k))],
  parameters = {l},
  variable = x,
  constraints = { x::positive },
  function = { functions:-cnormalI, functions:-erfcI, erfc },
  lhs = functions:-cnormalI(l,x),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "SM.cnormalI.sfrac.01",
  booklabel = "18.3.4",
  begin = [[1, x]],
  general = [[l+m-1, x]],
  parameters = {l},
  variable = x,
  constraints = { x::positive },
  function = { functions:-cnormalI, functions:-erfcI, erfc },
  lhs = functions:-cnormalI(l, x) / functions:-cnormalI(l-1, x),
  category = "S-fraction"
):

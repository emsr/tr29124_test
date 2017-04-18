create( 'series',
  label = "EX.ei.power.04",
  booklabel = "14.2.14",
  front = gamma + ln( x ),
  general = [ x^k / ( k * k! ) ],
  variable = x,
  constraints = { x::positive },
  function = Ei,
  lhs = Ei(x),
  category = "power series"
):

create( 'series',
  label = "EX.ein.power.01",
  booklabel = "14.2.16",
  factor = -1, 
  general = [ (-1)^(k+1) * z^(k+1) / ( (k+1) * (k+1)! ) ],
  function = { functions:-Ein, Ei },
  lhs = functions:-Ein(z),
  category = "power series"
):

create( 'series',
  label = "EX.ei.asymp.02",
  booklabel = "14.2.19",
  factor = exp(x) / x,
  general = [k! * x^(-k)],
  variable = x,
  function = Ei,
  lhs = Ei(x),
  category = "asymptotic series"
):

create( 'contfrac',
  label = "EX.ei.sfrac.02",
  booklabel = "14.2.21",
  front = 2 * sum(z^( 2*i + 1 ) / ( ( 2* i + 1 ) * ( 2 * i + 1 )! ), i=0..infinity),
  factor = -exp( -z ),
  begin = [[1 / z, 1]],
  general = [[ floor(m/2) / z, 1 ]],
  constraints = { x::positive },
  function = Ei,
  lhs = Ei(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EX.ein.sfrac.01",
  booklabel = "14.2.23",
  front = gamma + ln(z),
  factor = exp( -z ),
  begin = [[1/z, 1]],
  general = [[ floor(m/2) / z, 1]],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = { functions:-Ein, Ei },
  lhs = functions:-Ein(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EX.ei.cfrac.03",
  booklabel = "14.2.24",
  begin = [[exp(x) / x, 1]],
  general = [[ -floor(m/2) / x, 1]],
  variable = x,
  constraints = { x::positive },
  function = Ei,
  lhs = Ei(x),
  category = "C-fraction"
):

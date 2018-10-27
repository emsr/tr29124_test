create( 'series',
  label = "SM.regincbeta.power.01",
  booklabel = "18.5.16a",
  factor = x^a / (a * Beta(a,b)),
  general = [ pochhammer(a, k) * pochhammer(1-b, k) * x^k / (pochhammer(a+1, k) * k!) ],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "power series"
):

create( 'series',
  label = "SM.regincbeta.power.02",
  booklabel = "18.5.16b",
  factor = x^a * (1-x)^b / (a * Beta(a,b)),
  general = [ pochhammer(1, k) * pochhammer(a+b, k) * x^k / (pochhammer(a+1, k) * k!) ],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "power series"
):

create( 'contfrac',
  label = "SM.regincbeta.cfrac.01",
  booklabel = "18.5.17",
  factor = x^(a-1) * (1-x)^b / (a * Beta(a,b)),
  begin = [[x, 1]],
  even = [ - (a+(m-2)/2)*(a+b+(m-2)/2)/((a+m-2)*(a+m-1)) * x, 1],
  odd = [ ((m-1)/2)*(b-(m-1)/2)/((a+m-2)*(a+m-1)) * x, 1],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "C-fraction"
):

create( 'contfrac',
  label = "SM.regincbeta.mfrac.01",
  booklabel = "18.5.18",
  factor = x^a * (1-x)^b / Beta(a, b),
  begin = [[1, a + (1-a-b)*x]],
  general = [[ (m-1)*(b-m+1)*x, a+m-1 - (a+b-m)*x ]],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "M-fraction"
):

create( 'contfrac',
  label = "SM.regincbeta.norlund.01",
  booklabel = "18.5.19",
  factor = x^a * (1-x)^b / Beta(a, b),
  begin = [[1, a - (a+b)*x]],
  general = [[ (a+b+m-2)*(m-1)*(x-x^2), a+m-1 - (a+b+2*m-2)*x ]],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1/2)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "Norlund fraction"
):

create( 'contfrac',
  label = "SM.regincbeta.frac.01",
  booklabel = "18.5.20",
  factor = x^a * (1-x)^b / Beta(a, b),
  begin = [[1, (a * (a+1-(a+b)*x)) / (a+1)]],
  general = [[(a+m-2)*(a+b+m-2)*(m-1)*(b-m+1)*x^2/(a+2*m-3)^2, (m-1)+(m-1)*(b-m+1)*x/(a+2*m-3)+(a+m-1)/(a+2*m-1) *(a+1-(a+b)*x+(m-1)*(2-x))]],
  parameters = {a, b},
  variable = x,
  constraints = { a::positive, b::positive, x::RealRange(0, Open(1)) },
  function = { functions:-regincBeta, functions:-incBeta, Beta },
  lhs = functions:-regincBeta(x, a, b),
  category = "continued fraction"
):

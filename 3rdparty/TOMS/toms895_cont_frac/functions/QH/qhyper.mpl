create( 'contfrac',
  label = "QH.qhyper.cfrac.01",
  booklabel = "19.2.1",
  front = 1,
  even = [ (1 - b*q^(m/2)) * (c*q^(m/2) - a) * q^(m/2 - 1) / ( (1 - c*q^(m - 1)) * (1 - c*q^m) ) * z, 1 ],
  odd = [ (1 - a*q^((m-1)/2)) * (c*q^((m-1)/2) - b) * q^((m-1)/2) / ( (1 - c*q^(m - 1)) * (1 - c*q^m) ) * z, 1 ],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a,b*q],[c*q],q,z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "QH.qhyper.cfrac.02",
  booklabel = "19.2.2",
  begin = [[ 1 , 1 ]],
  even = [ (1 - a*q^((m-2)/2)) * (c*q^((m-2)/2) - 1) * q^((m-2)/2) / ( (1 - c*q^(m - 2)) * (1 - c*q^(m-1)) ) * z, 1 ],
  odd = [ (1 - q^((m-1)/2)) * (c*q^((m-1)/2) - a) * q^((m-1)/2 - 1) / ( (1 - c*q^(m - 2)) * (1 - c*q^(m-1)) ) * z, 1 ],
  parameters = {a, c, q},
  constraints = { abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,q],[c*q],q,z),
  category = "C-fraction"
):

create( 'contfrac',
  label = "QH.qhyper.tfrac.01",
  booklabel = "19.2.5a",
  front = (q*(1-c) + (a-b*q)*z) / (q*(1-c)),
  factor = 1/(q*(1-c)),
  general = [[ q*(1-b*q^m)*(c*q^m-a) * z, q*(1-c*q^m) + (a-b*q^(m+1))*z ]],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))), abs(z) < abs(q/a) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a,b*q],[c*q],q,z),
  category = "T-fraction"
):

create( 'contfrac',
  label = "QH.qhyper.tfrac.02",
  booklabel = "19.2.5b",
  front = 1 + q*(1-c) / ( (a-b*q)*z ),
  factor = 1/( (a-b*q)*z ),
  general = [[ q*(1-b*q^m)*(c*q^m-a) * z, q*(1-c*q^m) + (a-b*q^(m+1))*z ]],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, (log[q](c)) :: Not(AndProp(integer, RealRange(-infinity, Open(0))) ), abs(z) > abs(q/a) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a,b*q],[c*q],q,z),
  category = "T-fraction"
):


create( 'contfrac',
  label = "QH.qhyper.norlund.01",
  booklabel = "19.2.7",
  front = (1-c - (a+b-a*b-a*b*q)*z) / (1-c),
  factor = 1/(1-c),
  general = [[ (1-a*q^m)*(1-b*q^m)*(c*z-a*b*q^m*z^2)*q^(m-1), (1-c*q^m) - (a+b-a*b*q^m-a*b*q^(m+1))*q^m*z ]],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))), abs(z) < abs(q/a) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a*q,b*q],[c*q],q,z),
  category = "Norlund fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.01",
  booklabel = "19.2.10",
  front = 1,
  general = [[a, -1], [a*q*(1-q^( (m+1)/3 - 1) * z), c-a*q^( (m+1)/3 )*z], [a*b*q^(m/3)*z-c,1] ],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))), abs(z) < 1 },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a*q,b],[c],q,z),
  category = "continued fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.02",
  booklabel = "19.2.11",
  begin = [[ (1-z), (1+c*q^(-1) -(a+b)*z) ]],
  general = [[ -(c*q^(-1)-a*b*q^(m-2)*z)*(1-q^(m-1)*z), 1+c*q^(-1)-(a+b)*q^(m-1)*z ]],
  parameters = {a, b, c, q},
  constraints = { abs(q) < 1, abs(c/q) <> 1 },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,q*z)/functions:-qhyper([a,b],[c],q,z),
  category = "continued fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.03",
  booklabel = "19.2.12",
  front = 1,
  general = [[ (1-a*q^( (m-1)/2 ) )*(1-b*q^( (m-1)/2) )* z, 1-z], [a*b*q^(m-1) * z - c * q^( m/2 - 1) , 1 ] ],
  parameters = {a, b, c, q},
  constraints = { z <> 1, abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b],[c],q,z)/functions:-qhyper([a,b],[c],q,q*z),
  category = "continued fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.04",
  booklabel = "19.2.13",
  front = 1,
  general = [[ (1-b*q^( (m-1)/2 ))* z, (1-a)*(1-z)], [a - c * q^( m/2 - 1) , 1 ] ],
  parameters = {a, b, c, q},
  constraints = { abs( z/((1-a)*(1-z)) ) < 1/4, abs(a/((1-a)*(1-z))) < 1/4, abs(q) < 1, (log[q](c))::Not(AndProp(integer,RealRange(-infinity,Open(0)))) },
  function = functions:-qhyper,
  lhs = functions:-qhyper([a*q,b],[c],q,z)/functions:-qhyper([a,b],[c],q,q*z),
  category = "continued fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.05",
  booklabel = "19.3.1",
  front = 1,
  even = [- (1-e*q^(((m/2)-1)+1)/a)*(1-e*q^(((m/2)-1)+1)/b)*(1-e*q^(((m/2)-1)+1)/c)*f*q^((m/2)-1) / ((1-e*q^(2*((m/2)-1)+1))*(1-e*q^(2*((m/2)-1)+2))*(1-f*q^((m/2)-1))), 1 ],
  odd = [ (e*f/(a*b*c)) * (1-a*q^((m-1)/2))*(1-b*q^((m-1)/2))*(1-c*q^((m-1)/2))/((1-e*q^(2*((m-1)/2)))*(1-e*q^(2*((m-1)/2)+1))*(1-f*q^((m-1)/2))) , (1-e*f/(a*b*c)) / (1- f*q^((m-1)/2)) ],
  parameters = {a, b, c, e, f, q},
  function = functions:-qhyper,
  lhs = functions:-qhyper([a,b,c],[e,f],q,[e*f/(a*b*c)])/functions:-qhyper([a,b,c],[e*q,f],q,[e*f*q/(a*b*c)]),
  category = "continued fraction"
):

create( 'contfrac',
  label = "QH.qhyper.frac.06",
  booklabel = "19.3.2",
  front = 1,
  even = [ a * (1 - e*q^( (m/2) - 1 ) / a ) * ( 1 - f*q^((m/2)-1) / a), 1 ],
  odd = [ (1- b*q^((m-1)/2)) * (1-c*q^((m-1)/2)) * e*f / (a*b*c*q), (1-a)*(1-e*f/(a*b*c*q)) ],
  parameters = {a, b, c, e, f, q},
  function = functions:-qhyper,
  lhs = functions:-qhyper([a*q,b,c],[e,f],q,[e*f/(a*b*c*q)])/functions:-qhyper([a,b,c],[e,f],q,[e*f/(a*b*c)]),
  category = "continued fraction"
):

create( 'series',
  label = "EF.tan.power.01",
  booklabel = "11.3.3",
  general = [ (4^(k+1)*(4^(k+1)-1)*abs(bernoulli(2*(k+1)))/(2*(k+1))!)*z^(2*(k+1)-1) ],
  constraints = { abs(z) < Pi/2 },
  function = tan,
  lhs = tan(z),
  category = "power series"
):

create( 'contfrac',
  label = "EF.tan.sfrac.01",
  booklabel = "11.3.7",
  begin = [[z, 1]],
  general = [ [-z^2/((2*m-1)*(2*m-3)), 1] ],
  function = tan,
  lhs = tan(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.tan.thiele.01",
  booklabel = "11.3.8",
  begin = [[z, 1]],
  general = [ [(2*m-3)^2-z^2, 2] ],
  function = tan,
  lhs = tan(Pi*z/4),
  category = "Thiele interpolating fraction" 
):

create( 'contfrac',
  label = "EF.tan.thiele.02",
  booklabel = "11.3.9",
  begin = [[z, 1], [-4*Pi^(-2)*z^2,1]],
  general = [ [(m-2)^4-4*Pi^(-2)*(m-2)^2*z^2, 2*m-3] ],
  function = tan,
  lhs = tan(z),
  category = "Thiele interpolating fraction"
):

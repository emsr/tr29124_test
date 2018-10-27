create( 'series',
  label = "GA.binet.asymp.01",
  booklabel = "12.2.6",
  general = [ ( ( bernoulli( 2 * ( k + 1 ) ) ) / ( ( 2 * k + 1 ) * ( 2 * k + 2) ) ) * z^( -2 * k - 1 )  ],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = functions:-Binet,
  lhs = functions:-Binet(z),
  category = "asymptotic series"
):

$ifdef QD
create( 'contfrac',
  label = "GA.binet.sfrac.01",
  booklabel = "12.2.11",
  general = [ [qdbinet(m), z] ],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = functions:-Binet,
  lhs = functions:-Binet(z),
  category = "S-fraction"
):
$endif

$ifdef QD
create( 'contfrac',
  label = "GA.binet.sfrac.02",
  booklabel = "12.2.12",
  factor = 1 / sqrt( z ),
  general = [ [qdbinet(m) * z, 1] ],
  constraints = { abs(functions:-argument(z)) < Pi },
  function = functions:-Binet,
  lhs = functions:-Binet(1/sqrt(z)),
  category = "S-fraction"
):
$endif

create( 'series',
  label = "CN.rabbit.power.01",
  booklabel = "10.10.2",
  general = [ 2^( -floor( ( k + 1 )* ( (sqrt(5) + 1 ) / 2 ) ) ) ],
  function = functions:-rabbit,
  lhs = functions:-rabbit(),
  category = "power series"
):

create( 'contfrac',
  label = "CN.rabbit.rfrac.01",
  booklabel = "10.10.5",
  general = [ [1, 2^( :-combinat[fibonacci](m - 1) ) ] ],
  function = functions:-rabbit,
  lhs = functions:-rabbit(),
  category = "regular continued fraction"
):

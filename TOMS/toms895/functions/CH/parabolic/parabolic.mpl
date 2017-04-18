create( 'contfrac',
  label = "CH.cylinderU.cfrac.01",
  booklabel = "16.5.7",
  begin = [[1, x]],
  general = [ [ a+(m-1)-1/2, x]],
  parameters = {a},
  variable = x,
  constraints = { x::positive, (2*a)::posint },
  function = CylinderU,
  lhs = CylinderU(a, x) / CylinderU(a-1, x),
  category = "C-fraction"
):

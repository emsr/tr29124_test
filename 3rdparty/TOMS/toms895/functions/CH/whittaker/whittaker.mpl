create( 'series',
  label = "CH.whittakerW.asymp.01",
  booklabel = "16.4.7",
  factor = exp(-z^2) * z^kappa,
  general = [pochhammer(-kappa - mu + 1/2, k) * pochhammer(-kappa + mu + 1/2, k) * (-z)^(-k) / k!],
  parameters = {kappa, mu},
  constraints = { abs(functions:-argument(z)) < 3*Pi/2 },
  function = WhittakerW,
  lhs = WhittakerW(kappa, mu, z),
  category = "asymptotic series"
);

create( 'series',
  label = "CH.whittakerPsi.asymp.01",
  booklabel = "16.4.12",
  general = [pochhammer(alpha + 1/2, k) * (beta + 1/2, k) * z^(-k-1) / k! ],
  parameters = {alpha, beta},
  function = functions:-WhittakerPsi,
  lhs = functions:-WhittakerPsi(alpha, beta, z),
  category = "asymptotic series"
);


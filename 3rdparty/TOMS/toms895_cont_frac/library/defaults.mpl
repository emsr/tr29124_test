##### ============================
##### ===== DEFAULT SETTINGS =====
##### ============================

# this table defines the default values and colors that
# will be used when they are not specified explicitly

defaults := table([

  terms_series = 5,
  terms_contfrac = 5,

  table_terms_series = 10,
  table_terms_contfrac = 10,

  plot_color_function = blue,
  plot_color_approx = green,
  plot_color_error = red,

  # if you want to make it possible for users to add formulas
  # to the formula table, then you have to set the variable
  # 'defaults[add_formulas]' to 'true' ...

  add_formulas = false,

  # the following variable defines the maximum amount of steps
  # that are performed after which a simplification is done;
  
  evaluate_simplify = 2

]);

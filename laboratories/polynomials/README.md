
These are experiments with NRiCs rational polynomial function fitter.
They, and others, just do a linear fir to P(x) - f(x)Q(x) = 0;

* [ ]  Make extrema a product of the linear fit since that routine oversamples the function anyway.
* [ ]  Get a Remes to take the linear fit output and generate a real minimax rational.
* [ ]  Make a switch that lets you weight by function value to get relative error rather than absolute error.
* [ ]  Make a switch to set q_n to one (monic denominator polynomial) rather than q_0 = 1.

I think having a linear fit as an accessible tool is just fine.
Then optionally feed that output to Remes.

Having Remes search for minimax points using optimization tools would be scary.

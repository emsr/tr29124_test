# TODO List for cxx_differentiation

* [ ] Add Doxygen and comment everything nicely.

* [ ] Add initial stepsize heuristics.
      Idea - call central and forward thins *_step.

* [ ] See if these routines could be made to work with complex or vector types.
      The struct will need separate types for the absolute-value error estimates.
      This is like a gradient component thing.

* [ ] Make demos that shows all of the above plus derivatives of functions
      of a complex or vector returning a real scalar

* [ ] WWBD? What would Boost do?  Look at their APIs and functionality.
      Apparently nothing yet.

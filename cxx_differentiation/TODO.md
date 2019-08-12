# TODO List for cxx_differentiation

*   Add Doxygen and comment everything nicely.

*   Add initial stepsize heuristics.
    Idea - call central and forward thins *_step.
    This was tried and found to be of dubious use.
    Mybe try again with bad stepsize choices.

*   See if these routines could be made to work with complex or vector types.
    The struct will need separate types for the absolute-value error estimates.
    This is like a gradient component thing.
    I need a tool to extract the scalar type either by detection of value_type
    or decltype(abs(x)) or decay_t.

*   Make demos that shows all of the above plus derivatives of function
    of a complex or vector returning a real scalar

*   WWBD? What would Boost do?  Look at their APIs and functionality.
    Apparently nothing yet.
